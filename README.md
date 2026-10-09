# AutoATES-v2.0
A fully automated algorithm used to provide ATES mapping for large areas.

## Overview

`AutoATES_classifier.py` combines five raster layers into an Avalanche Terrain
Exposure Scale (ATES) map with classes 0 to 4:

1. Slope angle, derived from the DEM.
2. Runout zones, from the Flow-Py travel-angle layer.
3. Avalanche frequency proxy, from the Flow-Py cell-count layer.
4. Forest cover, from a canopy layer.
5. Potential release areas (PRA), from a binary start-zone mask.

The combined classification is then generalized by removing clusters smaller
than a configured area.

TODO: describe how this repository relates to the upstream project at
<https://github.com/AutoATES/AutoATES-v2.0> and what, if anything, differs.

## Requirements

- conda (dependencies come from `conda-forge`)
- Python 3.9 or newer
- numpy, scipy, rasterio, gdal, scikit-image, pyyaml

The full list is in `environment.yml`.

## Setup

```bash
conda env create -f environment.yml
conda activate autoates_env
```

## Usage

1. Place the input rasters in the working directory (see `data/README.md`).
2. Review `config.yaml`.
3. Run the classifier:

```bash
python AutoATES_classifier.py
```

By default the script reads `config.yaml` from its own folder. To use a
different file:

```bash
python AutoATES_classifier.py --config path/to/config.yaml
```

Outputs are written to the working directory. The final map is `ates_gen.tif`.

TODO: document how the input layers (PRA, Flow-Py travel angle and cell count)
are produced; that preprocessing is not part of this repository.

## Configuration

All parameters live in `config.yaml`:

| Section | Purpose |
|---|---|
| `working_dir` | Folder holding the inputs and receiving the outputs. Relative paths resolve against the config file's folder. |
| `inputs` | Input raster file names, resolved against `working_dir`. |
| `logging` | Log level and log file name. |
| `nodata_value` | NoData value written to every output raster. |
| `forest_type` | Canopy layer type: `bav`, `stems`, `pcc` or `sen2cc`. Selects the matching `tree_thresholds` entry. |
| `slope_angle_thresholds`, `smoothing_window_size` | Slope class boundaries and the smoothing window used for the class 4 threshold. |
| `alpha_angle_thresholds` | Travel-angle boundaries for runout classes 1 to 3. |
| `tree_thresholds` | Forest density boundaries per `forest_type`. |
| `cell_count_thresholds` | Cell-count class boundaries. |
| `class_codes`, `ates_lookup` | Lookup of the final ATES class for each release / forest / terrain combination. |
| `generalization` | Minimum cluster size and fill settings for the final map. |

## Structure

```
.
├── AutoATES_classifier.py   classifier script
├── config.yaml              classifier parameters
├── environment.yml          conda environment
├── data/                    inputs and outputs (git-ignored, see data/README.md)
├── docs/                    further documentation (empty)
├── tests/                   tests (empty)
├── CITATION.cff
└── LICENSE
```

## Tests

No tests exist yet. `tests/` is a placeholder.

## Branching

`main` is the only long-lived branch. Other branches are short-lived, branch off
`main`, and are named `<type>/<short-kebab-description>`, where type is one of
`feat`, `fix`, `refactor`, `docs`, `data`, `exp` or `chore`.

## Citation

See `CITATION.cff`.

## License

GPL-3.0. See `LICENSE`.

## Maintainer

TODO: name the maintainer and a contact.
