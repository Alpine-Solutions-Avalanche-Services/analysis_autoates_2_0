# `data/` — AutoATES inputs and outputs

This folder holds the raster inputs consumed by `AutoATES_classifier.py` and the
intermediate and final layers it generates. **The contents of this folder are
git-ignored** (only this README is tracked) because the rasters are large and
case-specific. Place your own data here following the structure below.

## Required input layers

These are the layers the classifier reads. Their file names are set under
`inputs` in `config.yaml`:

| File | Variable | Description |
|------|----------|-------------|
| `dem.tif` | `DEM` | Digital elevation model. Drives slope-angle classification. |
| `forest_ccp.tif` | `canopy` | Forest canopy layer. Units depend on `forest_type` (`bav`, `stems`, `pcc`, `sen2cc`). |
| `cellCounts.tif` | `cell_count` | Flow-Py cell-count layer (avalanche frequency proxy). |
| `fpTravelAngle.tif` | `FP` | Flow-Py travel-angle (alpha angle) layer used for runout-zone classification. |
| `PRA_binary.tif` | `SZ` | Binary potential release area / start-zone mask (1 = release area, 0 = otherwise). |

## Generated layers

The classifier writes these into the working directory as it runs:

| File | Stage |
|------|-------|
| `slope.tif` | Reclassified slope-angle terrain classes (0–4). |
| `slope_smooth.tif` | Slope smoothed with the `WIN_SIZE` moving window (used for the Class 4 threshold). |
| `flowpy.tif` | Runout zones derived from the travel-angle thresholds (AAT1–AAT3). |
| `cellcount_reclass.tif` | Cell-count layer reclassified to classes 1–3 (CC1, CC2). |
| `forest_reclass.tif` | Canopy reclassified into open / sparse / dense / very dense (TREE1–TREE3). |
| `SZ_reclass.tif` | Start-zone mask rescaled for the combination step. |
| `merge_new.tif` | Maximum of slope, runout, and cell-count classes. |
| `merge_all.tif` | Combined terrain, forest, and start-zone classification before generalization. |
| `ates_gen.tif` | **Final output:** generalized ATES map (clusters smaller than `ISL_SIZE` removed). |
| `inputpara.csv` | Log of the input parameters used for the run (appended each run). |
| `autoates.log` | Run log (appended each run). The name is set under `logging` in `config.yaml`. |

## Subfolders

- `PRA/` — outputs of the potential-release-area / Flow-Py preprocessing step
  (`PRA_binary.tif`, `PRA_continous.tif`, `windshelter.tif`, `log.txt`).
- `raster_data/` — a self-contained example/working set mirroring the layers above.

## Notes

- GIS sidecar files (`.aux.xml`, `.ovr`, `.tfw`, `.xml`, `.vat.dbf`, `.vat.cpg`)
  accompany many rasters and are produced/consumed by GIS software.
- ArcGIS `*.sr.lock` files are transient locks and should not be committed.
- The classifier's working directory is set by `working_dir` in `config.yaml`;
  point it at this folder (or a subfolder) before running.
