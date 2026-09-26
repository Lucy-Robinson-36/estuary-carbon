# estuary-carbon
Estuary-Scale Assessment of Intertidal Organic Carbon: A Multi-Method Approach

A workflow for calculating and modelling blue carbon (organic carbon) stocks in tidal flat and saltmarsh environments, from raw core/geochemical data through to spatial predictions.

## Set Up Files

1. Create a folder called `Carbon_Map`. This folder will be your working directory.
2. Within `Carbon_Map`, create four further folders:
   - `predictors`
   - `samples`
   - `scripts`
   - `outputs`

### Populate the `predictors` folder

- A `.tif` file for each predictor
- A `.shp` file for each AOI (Area of Interest)
- A `.shp` file for each of the following:
  - Single polygon — tidal flat
  - Single polygon — saltmarsh
  - Multi-polygon — individual marshes

### Populate the `samples` folder

Add CSV files for:

- Geochemical results
- Coordinates
- Bulk density
- Marsh horizons

**Notes on sample data:**

- Sample names are presumed to take the format `Site_Transect_Core_depth1-depth2`
- Geochemical results are presumed to come from an Elementar SoliTOC. The code can be adapted for other outputs (e.g. TGA data) but will need customising first
- Coordinates should be in `EPSG:4326`
- Bulk density — units g/cm3
- Marsh horizons — thickness of the fibrous peat, humified peat, and transitional sediment layers, with the following columns:
  - `Site`
  - `Transect`
  - `Core`
  - `Fibrous_peat_thickness_cm`
  - `Humified_peat_thickness_cm`
  - `Transitional_thickness_cm`

### Populate the `scripts` folder

Download the following scripts from GitHub and move them into the `scripts` folder:

- `Core_calculations`
- `Data_preparation`
- `Full_Surface_TOC`
- `Full_Surface_OCl`
- `Marsh_Surface_TOC`
- `Marsh_Surface_OCl`
- `Marsh_Depth_TOC`
- `Marsh_Depth_OCl`
- `Model_selection`
- `Final_models`
- `Model_predictions`

### Set up the `outputs` folder

Within `outputs`, create four further folders (these will be populated automatically as the code runs):

- `core_results`
- `model_inputs`
- `initial_models`
- `final_models`

## Workflow

Run the scripts in the following order:

1. **Core_calculations**
2. **Data_preparation**
3. **Initial Models**
   1. Full_Surface_TOC
   2. Full_Surface_OCl
   3. Marsh_Surface_TOC
   4. Marsh_Surface_OCl
   5. Marsh_Depth_TOC
   6. Marsh_Depth_OCl
4. **Model_selection**
5. **Final_models**
6. **Model_predictions**

## Folder Structure

```
Carbon_Map/
├── predictors/
│   ├── *.tif                  # Predictor rasters
│   ├── AOI.shp                # AOI boundary(ies)
│   ├── tidal_flat.shp         # Single polygon tidal flat
│   ├── saltmarsh.shp          # Single polygon saltmarsh
│   └── individual_marshes.shp # Multi-polygon individual marshes
├── samples/
│   ├── geochemical_results.csv
│   ├── coordinates.csv
│   ├── bulk_density.csv
│   └── marsh_horizons.csv
├── scripts/
│   ├── Core_calculations
│   ├── Data_preparation
│   ├── Full_Surface_TOC
│   ├── Full_Surface_OCl
│   ├── Marsh_Surface_TOC
│   ├── Marsh_Surface_OCl
│   ├── Marsh_Depth_TOC
│   ├── Marsh_Depth_OCl
│   ├── Model_selection
│   ├── Final_models
│   └── Model_predictions
└── outputs/
    ├── core_results/
    ├── model_inputs/
    ├── initial_models/
    └── final_models/
```
