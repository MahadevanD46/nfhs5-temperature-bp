# Measurement-day ambient temperature, blood pressure and blood glucose in NFHS-5

Code for: Duraiswamy M. *Measurement-day ambient temperature, blood pressure and blood glucose in India's
National Family Health Survey (NFHS-5, 2019–21): a cross-sectional analysis.* (under review, *International
Journal of Epidemiology*).

Persistent archive: DOI to be added after the first Zenodo release

## What is here

```
code/
  A_original_pipeline/            main analysis (as in the submitted manuscript)
    ERA_download_fullset.R          ERA5-Land hourly 2-m temperature -> daily means (UTC+05:00 local days)
    download_normals.R              ERA5-Land monthly means 1991-2018 (climatological normals)
    nfhs5_temperature_full_period.R NFHS-5 preparation, exposure linkage, main models
    prepare_district_boundaries.R   NFHS-5 district boundaries (DHS Spatial Data Repository), one polygon per district code
  B_revision/                     analyses added in response to peer review (run after A)
    00_config.R                     paths, model specification, helpers  <- edit paths here only
    01_flow_validation.R            participant flow, excluded vs included, BP validation, date handling
    02a_download_humidity_rain.R    ERA5-Land dewpoint + precipitation (needs CDS key)
    02b_extra_exposures.R           interview-date and displacement-weighted temperature, humidity,
                                    coastal-infill method, grid cell -> analysis_dataset_v2.rds
    02c_timezone_check.R            hourly vs daily, UTC+05:00 vs IST check
    02d_imd_validation.R            OPTIONAL: comparison with IMD gridded station temperature
    03_descriptives_diagnostics.R   Table 1, state x month counts, date gaps, shifts, identifying variation
    04_models.R                     sensitivity analyses (all outcomes), Conley/two-way SEs, splines,
                                    interaction tests, placebo diagnostics, state heterogeneity
    04c_additions.R                 counts in the model sample, glucose denominator, minimally adjusted models
    05_standardisation.R            model-based standardisation + 8 alternative specifications
    06_bootstrap.R                  Rao-Wu design-based bootstrap (slopes re-estimated each replicate)
    07_lisa.R                       district LISA + bootstrap reproducibility, maps
    08_figures_basic.R              basic plots of the main results (default styling)
    08b_maps_basic.R                plain maps of the state results and district LISA classes
    run_all.R                       runs 01-08 in order, each in a fresh R session
docs/                             data dictionary and data-sharing statement
tests/make_synthetic_data.R       synthetic NFHS-5-like dataset (same file/variable names) for testing
```

## Data (not included)

| Source | Files | How to obtain |
|---|---|---|
| NFHS-5 (India DHS 2019-21) | `IAPR7EFL.DTA`, `IAGE7AFL.shp` | Register at https://dhsprogram.com and request the India 2019-21 survey including GPS data |
| District boundaries | `sdr_subnational_boundaries2.shp` (NFHS-5 districts) | https://spatialdata.dhsprogram.com |
| ERA5-Land | hourly / daily-statistics 2-m temperature, dewpoint, precipitation; monthly means 1991-2018 | Copernicus Climate Data Store (free account + API key); downloaded by the scripts |
| IMD gridded temperature (optional) | `Maxtemp_MaxT_YYYY.GRD`, `Mintemp_MinT_YYYY.GRD` | https://www.imdpune.gov.in |

DHS terms of use do not permit redistribution of the microdata, GPS coordinates or any person- or
cluster-level file derived from them. See `docs/DATA_SHARING_STATEMENT.md`.

## How to reproduce

1. Put the DHS files in the project folder (the folder that contains `code/`).
2. Set your Copernicus key **in the R session only** (never in a file):
   `Sys.setenv(CDS_KEY = "your-token")`.
3. Run the original pipeline: `code/A_original_pipeline/ERA_download_fullset.R`, `download_normals.R`,
   then `nfhs5_temperature_full_period.R` (creates `output_full/`) and `prepare_district_boundaries.R`
   (creates the district boundary file; set its path in `00_config.R`).
4. Run the revision analyses: `source("code/B_revision/run_all.R")` (creates `output_revision/`;
   models ~3 h, bootstrap ~3-4 h for 500 replicates).
5. Estimates are printed to the console (and to `output_revision/*.log` when run via `run_all.R`) and
   written as plain CSV files to `output_revision/`; basic plots go to `output_revision/figures/`.

Run with R 4.6.1 (fixest 0.14.2, data.table, terra, sf, spdep, metafor, Hmisc, ggplot2; exact versions
are reported in the article) and tested end to end on the synthetic dataset (`tests/make_synthetic_data.R`).

### Test without DHS data
```r
source("tests/make_synthetic_data.R")        # creates synthetic_project/
setwd("synthetic_project"); dir.create("code"); file.copy("../code", ".", recursive = TRUE)
source("code/A_original_pipeline/nfhs5_temperature_full_period.R")
source("code/B_revision/run_all.R")
```

## Licence
Code: MIT licence (see `LICENSE`). Data are subject to the DHS Program and Copernicus licences.
