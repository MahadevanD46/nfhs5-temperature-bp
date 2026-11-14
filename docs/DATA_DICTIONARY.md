# Data dictionary — analytic variables

Person-level analytic file: `output_full/analysis_dataset.rds` (original pipeline) and
`output_revision/analysis_dataset_v2.rds` (adds the revision variables). One row per non-pregnant adult
aged ≥18 years with valid blood pressure, a valid measurement date and a temperature value.
**Restricted — never share.**

## Identifiers and design

| Variable | Source | Definition |
|---|---|---|
| hv001 | PR | Cluster number (= DHSCLUST in GPS file) |
| hv002, hvidx | PR | Household number, line number (person key = hv001 + hv002 + hvidx) |
| hv021 | PR | Primary sampling unit |
| hv022 | PR | Sample stratum (district × urban/rural) |
| hv024 | PR | State/UT code (labels in `00_config.R`, `state_lab`) |
| shdist | PR | District code |
| region | derived | NFHS region (North, Central, East, Northeast, West, South) |
| shphase | PR | Survey phase (1, 2) |
| wt | derived | hv005 / 10^6 (household weight; not normalised, not recalibrated) |

## Dates

| Variable | Definition |
|---|---|
| bio_date | Biomarker visit date from hv807y (2019–2021), hv807m (1–12), hv807d (1–31); NA if invalid |
| int_date | Household interview date from hv007, hv006, hv016 |
| date | Measurement date used for exposure = bio_date, else int_date |
| date_fallback | 1 if bio_date missing/invalid and int_date used |
| gap | date − int_date (days) |
| gap_cat | |gap|: 0, 1, 2–3, 4–7, >7 days |
| same_day | 1 if bio_date = int_date |
| ym, dow | Year-month and day of week of `date` (fixed effects) |
| season | winter (Dec–Feb), pre-monsoon (Mar–May), monsoon (Jun–Sep), post-monsoon (Oct–Nov) |

## Outcomes

| Variable | Definition |
|---|---|
| sbp, dbp | Mean of valid readings 2 and 3 (shb25s/d, shb29s/d; single valid reading if only one); valid range SBP 70–270, DBP 40–150 mmHg; values on a 0.5-mmHg grid |
| htn_measured | Measured elevated BP: sbp ≥140 or dbp ≥90 (medication not counted) |
| htn_report | Composite: sbp ≥140 or dbp ≥90 or bp_med = 1 |
| bp_med | shb21 = 1 (current antihypertensive medication; skipped = 0) |
| glu | Random capillary glucose shb74, 40–500 mg/dL (995–998 = not measured → NA) |
| glu_med | shb57 = 1 (glucose-lowering medication) |
| height | ha3/hb3 ÷ 10, 100–220 cm (women 15–49, men 15–54) |
| bmi | ha40/hb40 ÷ 100, 12–60 |

## Covariates

| Variable | Definition |
|---|---|
| hv104 | Sex (1 men, 2 women) |
| hv105, age_grp | Age in years; groups (17,29], (29,39], … (69,Inf] = 18–29 … ≥70 |
| hv270 | Wealth quintile |
| hv025 | Residence (1 urban, 2 rural) |
| sh25, sh26 | Current tobacco use; alcohol use |
| shb14a, shb14b, shb14c | Ate; drank tea/coffee; smoked in the 30 min before BP measurement |
| ate_cat | Hours since last ate (shb53): <1, 1–2, 3–4, 5–8, >8 (glucose models) |
| bp_hour, hour_cat | Time of first BP reading (shb17, hhmm → hours, 6–21 valid); six bands + "missing" |
| hv018 | Interviewer code (interviewer FE within state: hv024^hv018) |
| fan, cooler | Household owns electric fan (sh50g) / air cooler (sh50q) |
| hv040 | Cluster altitude (m) |
| sub_bmi | In BMI/height subsample (women 15–49, men 15–54) |

## Exposures (ERA5-Land, 0.1°)

| Variable | Definition |
|---|---|
| t0 | Daily mean 2-m temperature (°C) on `date` at the cluster's cell (coastal clusters: 3×3 then 5×5 land-pixel mean); local day UTC+05:00 (= IST day for hourly data) |
| filled | Cluster temperature from coastal infilling |
| coast_method | own cell / 3x3 / 5x5 |
| t_lag1, t_lag2, t_mean02 | Temperature 1 and 2 days before; mean of days 0–2 |
| t_lead7 | Temperature 7 days after measurement (placebo exposure) |
| t_int | Temperature on the interview date |
| t_disp | Displacement-weighted temperature (average over cells weighted by simulated displacement probabilities) |
| p_same_cell, n_cells | Probability that the true location is in the assigned cell; number of cells reached by simulated locations |
| cell, cell_date | ERA5-Land cell id; cell × date (variance clustering) |
| lat, lon | Displaced DHS coordinates (restricted; used only for Conley SEs) |
| d2m, rh | Daily mean dewpoint (°C); relative humidity (%, Magnus formula) |
| tp_mm, tp_cat | Daily precipitation (mm); none / very light / light / moderate / heavy+ |
| t_annual | Cluster 1991–2018 annual mean (mean of 12 monthly normals) |
| t_doy | Calendar-day normal (two-harmonic fit) = t0 − anom0 |
| anom0 | Temperature anomaly = t0 − t_doy |
| t_season_norm | Mean normal over the season of measurement |
| shift | t0 − t_annual (standardisation shift) |
| tbin | Temperature category (<10, 10–15, 15–20, 20–25 [ref], 25–30, >30 °C) |
