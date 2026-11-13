# Examining 400 years of epidemics and pandemics to estimate future burden

Analysis code for the paper by Horacio Chacon-Torrico and Christopher J. L. Murray (Institute for
Health Metrics and Evaluation, University of Washington).

One command reproduces every number and pipeline-generated figure in the paper.

## Reproduce

```r
install.packages(c(
  "tidyverse", "data.table", "config", "evir", "mev", "MASS",
  "scales", "gridExtra", "patchwork", "ggrepel", "cowplot", "viridis",
  "yaml", "svglite", "ragg", "future", "future.apply"
))
```

```bash
Rscript src/run_all.R
```

Requires R ≥ 4.1.0 (tested on 4.5.2, macOS and Linux) and the input files described below.

Results are written to a timestamped `output/run_all_<date>/` directory: every figure
(PNG and SVG), every table and the `_variables.yml` the manuscript reads, plus
`execution_log.csv`, `session_info.txt` and `run_summary.md`. The run is seeded
(`seed: 12345` in `config.yml`) and reproduces the published numbers and figures.
The `session_info.txt` of the run behind the published values is included as
`session_info_published_run.txt`.

Installing the packages takes about 10 minutes on a typical desktop. A full run takes
about 30 minutes (28–32 minutes measured on an Apple-silicon laptop); `run_sens.R`
alone takes about 21 minutes. No non-standard hardware is needed; a multi-core CPU
shortens the bootstrap.

## Input data

`data/pandemic_events.csv` holds the 77 events the pipeline reads, taken from the event
inventory published with the paper (Supplementary Data 1). Deaths are in thousands; `pop`
is the world population at the time of the event.

The population and Global Burden of Disease inputs are not included in this repository.
To run the pipeline, place `pop_forecast.csv`, `pop_draws.csv` and
`daly_death_ratio_draws.rds` in `data/`; see the Data availability statement of the paper.

## What it runs

`src/run_all.R` runs the pipeline end to end. The main steps:

| Script | Produces |
|---|---|
| `run_PoT.R` | Main mixture-model analysis — Figures 1–3, Extended Data Figures 1 and 6, burden and scenario tables |
| `build_figure_5_exceedance.R` | Figure 5, the cumulative-deaths exceedance curve |
| `run_validation.R` | Out-of-sample validation — Figure 4 |
| `run_sens.R` | Sensitivity sweeps over cutoff, threshold and start year — Extended Data Figure 3 |
| `bulk_diagnostics.R` | Bulk-component diagnostics and model comparison |
| `tail_diagnostics.R` | Tail threshold diagnostics — Extended Data Figure 2 |
| `run_appendix_tables.R` | Temporal model-comparison and breakpoint tables — Extended Data Table 1 |
| `run_hiv_sensitivity.R` | Sensitivity analysis excluding HIV/AIDS |
| `run_yaml_export.R` | `_variables.yml`, the numbers quoted in the text |

Extended Data Figures 4 (pipeline schematic) and 5 (event-selection flow) are drawn by hand and
are not pipeline outputs.

## Configuration

All parameters are pre-set in `config.yml` at publication values. See [CONFIG.md](CONFIG.md) for the
reference.

## License

BSD 3-Clause. See [LICENSE](LICENSE).
