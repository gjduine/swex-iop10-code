# swex-iop10-code

Analysis and plotting scripts, and WRF/WPS namelists, for:

> Duine, G.-J., et al. (in review). Lee-Slope Jet, Hydraulic Jump and Rotor Characteristics from Observations and Large-Eddy Simulations during SWEX IOP 10. Journal of Geophysical Research: Atmospheres.

The paper analyzes Intensive Observational Period 10 (IOP 10, 12–13 May 2022) of the
Sundowner Winds Experiment (SWEX) over the Santa Ynez Mountains, Santa Barbara, California.
It uses Doppler lidars, radiosondes, and surface stations, together with WRF simulations at
1 km grid spacing and large-eddy simulations (LES) at 111 m grid spacing (with original and
1 km-smoothed terrain).

## Contents

### Figures

| Figure | Script | Language | Description |
|---|---|---|---|
| 1 | `Fig01_studyArea.ipynb` | Python | Study area, terrain, and instrument locations |
| 2a | `Fig02a_terrain_elevation_1km.m` | MATLAB | Terrain of the 1 km domain |
| 2b | `Fig02b_terrain_elevation_111m.m` | MATLAB | Terrain of the 111 m domain |
| 2c | `Fig02c_terrain_elevation_111m_smth1km.m` | MATLAB | Terrain of the 111 m domain, smoothed to 1 km |
| 2d | `Fig02d_terrain_elevation_zoom_111m_minus_111m_smth1km.m` | MATLAB | Terrain difference, 111 m minus 111 m smoothed |
| 3 | `fig03_plot_IOP10_lidar_VAD_stare_v1.m` | MATLAB | VAD wind profiles and vertical stares at HQ, SBBG, and MFS-2 |
| 4 | `fig04_RHI_0_v2.ipynb` | Python | SBBG lidar range–height indicator (RHI) scans |
| 5 | `fig05_RHI_0_v2.ipynb` | Python | SBBG lidar RHI scans |
| 6 | `fig06_wavelet_multi_height_UWOW.ipynb` | Python | Wavelet analysis of observed vertical velocity (UWOW lidar, MFS-2) at multiple heights |
| 7a | `fig07a_plot_IOP10_SBBG_v2.m` | MATLAB | SBBG lidar time–height wind speed |
| 7b | `fig07b_plot_IOP10_UWOW_v2.m` | MATLAB | UWOW lidar (MFS-2) time–height wind speed |
| 7c–d | `fig07cd_plot_IOP10_vertProf_1km_10min.m` | MATLAB | UWOW observations vs. WRF 1 km |
| 7e–f | `fig07ef_plot_IOP10_vertProf_111m.m` | MATLAB | UWOW observations vs. LES 111 m |
| 7g–h | `fig07gh_plot_IOP10_vertProf_111m_smth1km.m` | MATLAB | UWOW observations vs. LES 111 m with 1 km-smoothed terrain |
| 8 | `fig08_wavelet_multi_height_WRF.ipynb` | Python | Wavelet analysis of simulated vertical velocity (111 m, time-series output) |
| 9 | `fig09_combineThreeSims_v2.ipynb` | Python | North–south cross sections for the three simulations, with Froude-number flow transitions |
| 10 | `fig10_froude_comp_LES_meso_v10.m` | MATLAB | Froude number and near-surface variables: observations vs. 1 km and 111 m simulations |
| 11 | `fig11_plot2Dwinds_threeTimes_IOP10_v2.py` | Python | 10 m wind speed and vectors for the three simulations at three times |
| S5 | `figS05_timeseries_rotor_stations_IOP10.py` | Python | Surface station time series near the SBBG and MFS-2 lidars (requires `plot_stations_IOP10.py`) |

### Model configuration

| File | Description |
|---|---|
| `namelist.wps` | WPS configuration for the 1 km (d01) and 111 m (d02) domains used for the LES |
| `namelist.input_545_24_z107` | WRF configuration for the 9/3/1 km nested run (107 vertical levels) |
| `namelist.input_2Doms_ndown_111m` | WRF configuration for the 111 m LES, driven by the 1 km run through `ndown` |
| `terrain_smoothing_interpolation.ipynb` | Replaces the 111 m terrain (`HGT_M` in `geo_em.d02.nc`) with the 1 km terrain interpolated to the 111 m grid |

The 111 m simulation with smoothed terrain uses the same `namelist.input_2Doms_ndown_111m`,
with the `geo_em.d02.nc` modified by `terrain_smoothing_interpolation.ipynb`.

Output paths in the namelists point to the authors' directories on NCAR Derecho and are kept
for the record. Set `geog_data_path` in `namelist.wps` to your own copy of the WPS geographical
data.

## Data

The observations are available from the NSF NCAR Earth Observing Laboratory (EOL) SWEX archive:

- Doppler lidars (ISS): NSF NCAR/EOL Integrated Sounding System Team (2023). SWEX: ISS Wind Lidar Data Products. Version 1.0. https://doi.org/10.26023/Q28P-EEBS-0Y0E
- Radiosondes (Montecito): 
- UWOW lidar: *[DOI]*
- SBBG lidar: *[DOI]*


SCE and RAWS surface stations are available at various websites. For the paper, Synoptic Data was used

The WRF and LES output is available upon request from the corresponding author

## Requirements

**WRF/WPS:** WRF v4.6.0, WPS v4.2.2.

**Python** (3.x): numpy, scipy, pandas, xarray, netCDF4, matplotlib, cartopy, wrf-python, pycwt.

**MATLAB** (R20xx), with these File Exchange packages on the MATLAB path:

- [export_fig](https://github.com/altmany/export_fig)
- [bluewhitered](https://www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered)
- [distinguishable_colors](https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors)

## Usage

File and output paths are hard-coded for the authors' systems (the UCSB voyager server and NCAR
Derecho). 

## Citation
If you use this code, please cite the paper above and this archive:

> Duine, G.-J. (2026). swex-iop10-code (Version v1.0) [Software]. Zenodo. https://doi.org/*[DOI]*

## License

*[e.g. MIT — see `LICENSE`]*

## Contact

Gert-Jan Duine, Earth Research Institute, University of California, Santa Barbara
(duine@eri.ucsb.edu)
