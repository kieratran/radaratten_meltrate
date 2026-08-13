# Ice-Penetrating Radar Attenuation & Basal Melt Rate Calculation

MATLAB scripts for estimating ice-penetrating radar attenuation rates and converting them to basal melt rates using a 1.5D thermal model.

**Author:** Kiera Tran

## Workflow

Run in this order:

1. **`attenuation_calculation.m`** — Loads radargram data, detects the firn/ice and bed-echo layer boundaries, and computes depth-resolved attenuation rates (dB/km) using both linear and piecewise fitting.
2. **`melt_calculation.m`** — Simulates the basal melt rate vs. attenuation rate relationship from the 1.5D thermal model, then inverts radar-observed attenuation rates (output of step 1) into basal melt rate estimates.

## Requirements

- MATLAB (tested version not specified — update this line with your version)
- Custom function library in `~/Ross_projects/functions` (referenced via `addpath(genpath(...))` in `attenuation_calculation.m`), including at minimum:
  - `geo_correction`
  - `piecewise_fit`
  - `smoothn`
  - `atten_calc`
  - `depth_ave`
  - `piecewise_unc`
  - `slopeSE`
  - `plot_radargram`
  - `atten2melt`
- Radar data in `CSARP_standard` format (e.g. from [Operation IceBridge](https://data.cresis.ku.edu/data/rds/))

## Data inputs

### `attenuation_calculation.m`
Loads a `.mat` file (`radar_data.mat` by default) containing:

| Variable | Description |
|---|---|
| `Data` | Radar power/amplitude matrix |
| `surfLayer`, `bedLayer` | Indices of surface and bed reflectors |
| `time` | Fast-time sample vector |
| `lat`, `lon` | Along-track coordinates |
| `plane_elevation` | Aircraft/platform elevation |
| `surfTime` | Surface returned time |

### `melt_calculation.m`
Loads `Model_inputs.mat` containing:

| Variable | Description |
|---|---|
| `thickness` | Ice thickness (m) |
| `velocity` | Ice velocity (m/yr), e.g. from MEaSUREs v2 |
| `dTdl` | Horizontal temperature gradient (K/yr) |
| `smb` | Surface mass balance (m/yr) |
| `Ts` | Surface temperature (K) |
| `basal_slope`, `surface_slope` | Slope effects (m/yr) |

It then loads `Attenuation` and `Attenuation_unc` (the output of `attenuation_calculation.m`, typically saved into `radar_data.mat`) to compute melt rates.

## Output

- `attenuation_calculation.m` produces `Attenuation` and `Attenuation_unc` (mean depth-resolved attenuation rate and uncertainty per trace).
- `melt_calculation.m` produces a `MeltRate` struct with mean (`.mu`) and uncertainty (`.unc`) fields for both Siple Dome (`.SD`) and Taylor Dome (`.TD`) ice chemistry assumptions.

## Notes

- `melt_calculation.m` warns that `tmd.Na_SD` / `tmd.Na_TD` can grow very large in memory — watch usage on large datasets.
- Segment lengths in the attenuation piecewise fit are swept from 20–50% of ice thickness (5% steps), plus a full-thickness (100%) linear fit.
