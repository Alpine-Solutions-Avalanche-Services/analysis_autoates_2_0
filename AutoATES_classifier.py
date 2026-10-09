import numpy as np
import rasterio, rasterio.mask
from osgeo import gdal
import argparse
import logging
import os
import sys
from typing import Any
from skimage import morphology
import csv
import scipy.ndimage
import yaml
from rasterio.fill import fillnodata


def parse_args() -> argparse.Namespace:
    """Parse the command-line arguments."""
    default_config = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.yaml')
    parser = argparse.ArgumentParser(description='AutoATES v2.0 classifier')
    parser.add_argument('--config', default=default_config, help='Path to the YAML configuration file')
    return parser.parse_args()


def configure_console_logging() -> None:
    """Send log messages to the console."""
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    root.addHandler(handler)


def add_file_logging(log_file: str, level: str) -> None:
    """Send log messages to a file as well and apply the configured level."""
    try:
        handler = logging.FileHandler(log_file, encoding='utf-8')
    except OSError as error:
        raise RuntimeError(f"Could not open log file '{log_file}': {error}") from error
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    root = logging.getLogger()
    root.setLevel(level)
    root.addHandler(handler)


def load_config(config_path: str) -> dict[str, Any]:
    """Read the YAML configuration file."""
    try:
        with open(config_path, 'r', encoding='utf-8') as stream:
            config = yaml.safe_load(stream)
    except OSError as error:
        raise RuntimeError(f"Could not read config file '{config_path}': {error}") from error
    except yaml.YAMLError as error:
        raise RuntimeError(f"Config file '{config_path}' is not valid YAML: {error}") from error
    if not isinstance(config, dict):
        raise RuntimeError(f"Config file '{config_path}' must contain a mapping of settings.")
    return config


def validate_config(config: dict[str, Any]) -> None:
    """Check that the configuration holds every setting the classifier needs."""
    required_keys = (
        'working_dir',
        'inputs.dem', 'inputs.canopy', 'inputs.cell_count', 'inputs.flow_path', 'inputs.start_zone',
        'logging.level', 'logging.file',
        'nodata_value', 'forest_type', 'smoothing_window_size',
        'slope_angle_thresholds.sat01', 'slope_angle_thresholds.sat12',
        'slope_angle_thresholds.sat23', 'slope_angle_thresholds.sat34', 'slope_angle_thresholds.max',
        'alpha_angle_thresholds.aat1', 'alpha_angle_thresholds.aat2',
        'alpha_angle_thresholds.aat3', 'alpha_angle_thresholds.max',
        'tree_thresholds',
        'cell_count_thresholds.cc1', 'cell_count_thresholds.cc2', 'cell_count_thresholds.max',
        'class_codes.forest.open', 'class_codes.forest.sparse',
        'class_codes.forest.dense', 'class_codes.forest.very_dense',
        'class_codes.forest_excluded', 'class_codes.non_release', 'class_codes.release',
        'ates_lookup.non_release', 'ates_lookup.release',
        'generalization.island_size', 'generalization.connectivity',
        'generalization.search_distance_divisor', 'generalization.smoothing_iterations',
    )
    for dotted_key in required_keys:
        node = config
        for part in dotted_key.split('.'):
            if not isinstance(node, dict) or part not in node:
                raise ValueError(f"Missing configuration key: '{dotted_key}'")
            node = node[part]

    forest_type = config['forest_type']
    tree_thresholds = config['tree_thresholds']
    if forest_type not in tree_thresholds:
        raise ValueError(
            f"forest_type '{forest_type}' has no entry in tree_thresholds "
            f"(available: {', '.join(tree_thresholds)})"
        )
    for tree_key in ('tree1', 'tree2', 'tree3'):
        if tree_key not in tree_thresholds[forest_type]:
            raise ValueError(f"Missing configuration key: 'tree_thresholds.{forest_type}.{tree_key}'")

    # Terrain classes 0-4, one lookup entry each
    terrain_class_count = 5
    for release_state in ('non_release', 'release'):
        for forest_class in config['class_codes']['forest']:
            row = config['ates_lookup'][release_state].get(forest_class)
            if not isinstance(row, list) or len(row) != terrain_class_count:
                raise ValueError(
                    f"ates_lookup.{release_state}.{forest_class} must list "
                    f"{terrain_class_count} ATES classes"
                )


def resolve_working_dir(working_dir: str, config_path: str) -> str:
    """Return the working directory as an absolute path, relative to the config file if needed."""
    config_dir = os.path.dirname(os.path.abspath(config_path))
    wd = os.path.normpath(os.path.join(config_dir, working_dir))
    if not os.path.isdir(wd):
        raise FileNotFoundError(f"Working directory '{wd}' does not exist (working_dir in '{config_path}')")
    return wd


def AutoATES(config: dict[str, Any], wd: str) -> None:
    """Classify avalanche terrain into ATES classes and write the rasters to the working directory."""
    logger = logging.getLogger(__name__)

    # --- Set Input Files
    inputs = config['inputs']
    DEM = os.path.join(wd, inputs['dem'])
    canopy = os.path.join(wd, inputs['canopy'])
    cell_count = os.path.join(wd, inputs['cell_count']) # replace with z_delta in next iteration
    FP = os.path.join(wd, inputs['flow_path'])
    SZ = os.path.join(wd, inputs['start_zone'])

    NODATA = config['nodata_value']

    # Moving window size to smooth slope angle layer for calcuation of Class 4 extreme
    WIN_SIZE = config['smoothing_window_size']

    # --- Define slope angle Thresholds
    # SAT34 is calculated on a smoothed raster layer, so the slope angle value is not representative of real world values
    slope_thresholds = config['slope_angle_thresholds']
    SAT01 = slope_thresholds['sat01']
    SAT12 = slope_thresholds['sat12']
    SAT23 = slope_thresholds['sat23']
    SAT34 = slope_thresholds['sat34']
    SLOPE_MAX = slope_thresholds['max']

    # --- Define alpha angle thresholds
    alpha_thresholds = config['alpha_angle_thresholds']
    AAT1 = alpha_thresholds['aat1']
    AAT2 = alpha_thresholds['aat2']
    AAT3 = alpha_thresholds['aat3']
    ALPHA_MAX = alpha_thresholds['max']

    # --- Add tree coverage criteria (upper bounds of "open", "sparse" and "mixed")
    tree_thresholds = config['tree_thresholds'][config['forest_type']]
    TREE1 = tree_thresholds['tree1']
    TREE2 = tree_thresholds['tree2']
    TREE3 = tree_thresholds['tree3']

    # --- Add cell count criteria
    cell_count_thresholds = config['cell_count_thresholds']
    CC1 = cell_count_thresholds['cc1']
    CC2 = cell_count_thresholds['cc2']
    CC_MAX = cell_count_thresholds['max']

    class_codes = config['class_codes']
    forest_codes = class_codes['forest']
    FOREST_EXCLUDED = class_codes['forest_excluded']
    release_codes = {'non_release': class_codes['non_release'], 'release': class_codes['release']}
    ates_lookup = config['ates_lookup']

    # --- Threshold for number of cells in a cluster to be removed (generalization)
    generalization = config['generalization']
    ISL_SIZE = generalization['island_size']

    # --- Write input parameters to CSV file
    labels = ['DEM', 'canopy', 'cell_count', 'FP', 'SAT01', 'SAT12', 'SAT23', 'SAT34', 'AAT1', 'AAT2', 'AAT3', 'TREE1', 'TREE2', 'TREE3', 'CC1', 'CC2', 'ISL_SIZE', 'WIN_SIZE']
    csvRow = [inputs['dem'], inputs['canopy'], inputs['cell_count'], inputs['flow_path'], SAT01, SAT12, SAT23, SAT34, AAT1, AAT2, AAT3, TREE1, TREE2, TREE3, CC1, CC2, ISL_SIZE, WIN_SIZE]
    csvfile = os.path.join(wd, "inputpara.csv")
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(labels)
        wr.writerow(csvRow)
    logger.info("Input parameters appended to %s", csvfile)

    # --- Calculate slope angle
    def calculate_slope(DEM):
        gdal.DEMProcessing(os.path.join(wd, 'slope.tif'), DEM, 'slope')
        with rasterio.open(os.path.join(wd, 'slope.tif')) as src:
            slope = src.read()
            profile = src.profile
        return slope, profile

    logger.info("Calculating slope classes from %s", DEM)
    slope, profile = calculate_slope(DEM)
    slope = slope.astype('int16')

    slope_nd = np.where(slope < 0, 0, slope)

    # Optional function to calculat class 4 slope using a neighborhood function - controlled by WIN_SIZE input parameter
    # If WIN_SIZE is set to 1 this function does not do anything to the SAT34 threshold calculation
    slope_smooth = scipy.ndimage.uniform_filter(slope_nd, size = WIN_SIZE, mode = 'nearest')

    # Update metadata
    profile.update({"driver": "GTiff", "nodata": NODATA, 'dtype': 'int16'})

    # Reclassify
    slope[np.where((0 < slope) & (slope <= SAT01))] = 0
    slope[np.where((SAT01 < slope) & (slope <= SAT12))] = 1
    slope[np.where((SAT12 < slope) & (slope <= SAT23))] = 2
    slope[np.where((SAT23 < slope) & (slope <= SLOPE_MAX))] = 3
    slope[np.where((SAT34 < slope_smooth) & (slope_smooth <= SLOPE_MAX))] = 4

    with rasterio.open(os.path.join(wd, "slope.tif"), 'w', **profile) as dst:
        dst.write(slope)

    with rasterio.open(os.path.join(wd, "slope_smooth.tif"), 'w', **profile) as dst:
        dst.write(slope_smooth)

    # --- Read Flow Path Data Once
    logger.info("Classifying runout zones from %s", FP)
    with rasterio.open(FP) as src:
        array = src.read(1)  # Read the first band of the raster
        profile = src.profile  # Extract metadata/profile information
        array = array.astype('int16')  # Convert the array to 16-bit integer

    # Update metadata
    profile.update({"driver": "GTiff", "nodata": NODATA, 'dtype': 'int16'})




    # --- AAT1 (Class 1 Runout Zone)
    flowpy1 = np.copy(array)
    flowpy1[np.where(flowpy1 < AAT1)] = 0  # Assign 0 to areas below `AAT1`
    flowpy1[np.where((flowpy1 >= AAT1) & (flowpy1 < ALPHA_MAX))] = 1  # Assign value 1 to areas meeting Class 1 criteria

    # --- AAT2 (Class 2 Runout Zone)
    flowpy2 = np.copy(array)
    flowpy2[np.where(flowpy2 < AAT2)] = 0  # Assign 0 to areas below `AAT2`
    flowpy2[np.where((flowpy2 >= AAT2) & (flowpy2 < ALPHA_MAX))] = 2  # Assign value 2 to areas meeting Class 2 criteria

    # --- AAT3 (Class 3 Runout Zone)
    flowpy3 = np.copy(array)
    flowpy3[np.where(flowpy3 < AAT3)] = 0  # Assign 0 to areas below `AAT3`
    flowpy3[np.where((flowpy3 >= AAT3) & (flowpy3 < ALPHA_MAX))] = 3  # Assign value 3 to areas meeting Class 3 criteria

    # --- Combine Layers
    # Use np.maximum() to create a composite raster where higher classes take precedence
    flowpy = np.maximum(flowpy1, flowpy2)
    flowpy = np.maximum(flowpy, flowpy3)

    # --- Reshape for rasterio write compatibility
    flowpy = flowpy.reshape(1, flowpy.shape[0], flowpy.shape[1])

    # --- Update Metadata and Save Result
    # Update the profile with the correct data type and nodata value, if necessary
    profile.update({"driver": "GTiff", "nodata": NODATA, 'dtype': 'int16'})  # Set nodata value explicitly if needed

    # --- Save Raster to File
    with rasterio.open(os.path.join(wd, "flowpy.tif"), 'w', **profile) as dst:
        dst.write(flowpy)


    # --- Add cell count criteria

    # --- Reclassify cell count criteria
    logger.info("Reclassifying cell counts from %s", cell_count)
    with rasterio.open(cell_count) as src:
        array = src.read()
        array = array.astype('int16')
        profile = src.profile

        # Update metadata
        profile.update({"driver": "GTiff", "nodata": NODATA, 'dtype': 'int16'})

        # Reclassify
        array[np.where(array <= 0)] = 0
        array[np.where((0 < array) & (array <= CC1))] = 1
        array[np.where((CC1 < array) & (array <= CC2))] = 2
        array[np.where((CC2 < array) & (array <= CC_MAX))] = 3

    with rasterio.open(os.path.join(wd, "cellcount_reclass.tif"), 'w', **profile) as dst:
        dst.write(array)

    # --- Combine Tree coverage, slope class and cell count

    src1 = rasterio.open(os.path.join(wd, "slope.tif"))
    src1 = src1.read()

    src2 = rasterio.open(os.path.join(wd, "flowpy.tif"))
    src2 = src2.read()

    src3 = rasterio.open(os.path.join(wd, "cellcount_reclass.tif"))
    src3 = src3.read()

    ates = np.maximum(src1, src2)
    ates = np.maximum(ates, src3)

    with rasterio.open(os.path.join(wd, "merge_new.tif"), 'w', **profile) as dst:
        dst.write(ates)

    # --- Add tree coverage criteria

    src1 = rasterio.open(os.path.join(wd, "merge_new.tif"))
    src1 = src1.read()
    profile.update({"driver": "GTiff", "nodata": NODATA, 'dtype': 'int16'})


    # --- Reclassify using the forest criteria
    logger.info("Reclassifying forest cover from %s (forest_type: %s)", canopy, config['forest_type'])
    forest = rasterio.open(canopy).read()
    forest_open=forest
    forest_open[forest_open > TREE1] = FOREST_EXCLUDED
    forest_open[(forest_open >= 0) & (forest_open <= TREE1)] = forest_codes['open']

    forest = rasterio.open(canopy).read()
    forest_sparse=forest
    forest_sparse[forest_sparse > TREE2] = FOREST_EXCLUDED
    forest_sparse[forest <= TREE1] = FOREST_EXCLUDED
    forest_sparse[(forest > TREE1) & (forest <= TREE2)] = forest_codes['sparse']

    forest = rasterio.open(canopy).read()
    forest_dense=forest
    forest_dense[forest_dense > TREE3] = FOREST_EXCLUDED
    forest_dense[forest_dense <= TREE2] = FOREST_EXCLUDED
    forest_dense[(forest_dense > TREE2) & (forest_dense <= TREE3)] = forest_codes['dense']

    forest = rasterio.open(canopy).read()
    forest_vdense=forest
    forest_vdense[forest_vdense < TREE3] = FOREST_EXCLUDED
    forest_vdense[forest_vdense >= TREE3] = forest_codes['very_dense']

    src2=np.maximum(forest_open, forest_sparse)
    src2=np.maximum(src2, forest_dense)
    src2=np.maximum(src2, forest_vdense)

    with rasterio.open(os.path.join(wd, "forest_reclass.tif"), 'w', **profile) as dst:
        dst.write(src2)

    # --- Add PRA criteria
    logger.info("Adding release areas from %s", SZ)
    src3 = rasterio.open(SZ)
    src3 = src3.read()

    src3[np.where(0 == src3)] = release_codes['non_release']
    src3[np.where(1 == src3)] = release_codes['release']

    with rasterio.open(os.path.join(wd, "SZ_reclass.tif"), 'w', **profile) as dst:
        dst.write(src3)

    # --- Look up the ATES class for each release / forest / terrain combination
    summed = np.sum([src1, src2, src3], axis=0)
    array = np.copy(summed)

    for release_state, release_code in release_codes.items():
        for forest_class, forest_code in forest_codes.items():
            for terrain_class, ates_class in enumerate(ates_lookup[release_state][forest_class]):
                array[np.where(summed == release_code + forest_code + terrain_class)] = ates_class
    array[np.where(array < 0)] = 0

    array = array.astype('int16')

    # --- Save raster to path
    with rasterio.open(os.path.join(wd, "merge_all.tif"), "w", **profile) as dest:
        dest.write(array)

    # --- Remove clusters of raster cells smaller than ISL_SIZE
    logger.info("Generalizing clusters smaller than %s map units squared", ISL_SIZE)
    raster = gdal.Open(DEM)
    gt =raster.GetGeoTransform()
    pixelSizeX = gt[1]
    pixelSizeY =-gt[5]
    num_cells = np.around(ISL_SIZE / (pixelSizeX * pixelSizeY))
    logger.debug("Minimum cluster size: %s cells", num_cells)
    # --- Open file
    src1 = rasterio.open(os.path.join(wd, "merge_all.tif"))
    src1 = src1.read(1)

    # --- Change values to prepare for morphology and rasterio.fill
    src1 = src1 + 1
    src1 = src1.reshape(1, src1.shape[0], src1.shape[1])

    # --- Same as region group in arcmap. Each cluster gets a value between 1 and num_labels (number of clusters)
    # 20210430 JS changed connectivity to 2
    lab, num_labels = morphology.label(src1, connectivity=generalization['connectivity'], return_num=True)

    rg = np.arange(1, num_labels+1, 1)

    # --- Loop through all clusters and assign all clusters with less then ISL_SIZE to the value 0 (set null)
    for i in rg:
        occurrences = np.count_nonzero(lab == i)
        if occurrences < num_cells:
            lab[np.where(lab == i)] = 0

    # --- Save as dtype int16
    lab = lab.astype('int16')

    search_dist = num_cells / generalization['search_distance_divisor']
    #search_dist = num_cells

    # --- This algorithm will interpolate values for all designated nodata pixels (marked by zeros) (nibble)
    data = rasterio.fill.fillnodata(src1, lab, max_search_distance=search_dist, smoothing_iterations=generalization['smoothing_iterations'])

    # --- Change values back to standardized way of plotting ATES (0, 1, 2, 3 and 4)
    data = data - 1
    data[np.where(data < 0)] = NODATA
    data = data.astype('int16')
    profile.update({"driver": "GTiff", "nodata": NODATA, 'dtype': 'int16'})

    # --- Save raster to path
    with rasterio.open(os.path.join(wd, "ates_gen.tif"), "w", **profile) as dest:
        dest.write(data)
    logger.info("Final ATES map written to %s", os.path.join(wd, "ates_gen.tif"))


def main() -> None:
    """Run the classifier with the configuration given on the command line."""
    logger = logging.getLogger(__name__)
    args = parse_args()
    configure_console_logging()
    try:
        config = load_config(args.config)
        validate_config(config)
        wd = resolve_working_dir(config['working_dir'], args.config)
        add_file_logging(os.path.join(wd, config['logging']['file']), config['logging']['level'])
        logger.info("Running AutoATES with %s in %s", args.config, wd)
        AutoATES(config, wd)
    except Exception:
        logger.exception("AutoATES run failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
