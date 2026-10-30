"""
plot2Dwinds_threeTimes_IOP10_v2.py
3 rows × 3 columns panel figure of 10-m wind speed + vectors for IOP10.
  Rows    : 1 km  |  111 m smth-1km  |  111 m
  Columns : 2200 PDT May 12  |  0310 PDT May 13  |  0500 PDT May 13

One shared colorbar at the bottom.
Latitude tick labels only on left column; longitude only on bottom row.
Panel labels (a)–(i) in each panel.

v2: PBL top for the 1 km run also from the thv + 1.5 K method
    (as the 111 m runs), instead of WRF PBLH.
"""

import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import pandas as pd
from netCDF4 import Dataset
import wrf
from wrf import (getvar, to_np, latlon_coords, get_cartopy)
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# ============================================================
# ── USER SETTINGS ───────────────────────────────────────────
# ============================================================

baseFolder_1km = '/glade/campaign/univ/ucsb0023/duine/wrf/wrf460/sundowners/SWEX2022/IOP10/'
run_1km        = 'run_1km_545_24_z107_2way/'
fName_1km      = baseFolder_1km + run_1km + 'wrfout_d03_2022-05-11_18:00:00'

baseFolder_sm  = ('/glade/campaign/univ/ucsb0023/duine/wrf/wrf460/sundowners/SWEX2022/IOP10/'
                  'run_1km_545_24_z107_2way/ndown_111m_km_opt2_2Dom_smth1km/')
baseFolder_111 = ('/glade/campaign/univ/ucsb0023/duine/wrf/wrf460/sundowners/SWEX2022/IOP10/'
                  'run_1km_545_24_z107_2way/ndown_111m_km_opt2_2Dom/')

wrffi  = 'wrfout_d04_2022-05-'
tfiles = ['13_04', '13_05', '13_06', '13_07', '13_08', '13_09',
          '13_10', '13_11', '13_12', '13_13', '13_14', '13_15',
          '13_16', '13_17']

# ── Target times ─────────────────────────────────────────────
# Derived from the time-index mapping in plot2Dwinds_threeSims_downstreamOnly.ipynb
# time_map[(fi, t)] → idt_range starts at 204 (= '13_04' step 0 = 21:00 PDT May 12)
# Each file covers 1 hour in 6 × 10-min steps.
#   2200 PDT May 12 = 05:00 UTC May 13 → tfiles[1]='13_05' t=0  idt=210
#   0310 PDT May 13 = 10:10 UTC May 13 → tfiles[6]='13_10' t=1  idt=241
#   0500 PDT May 13 = 12:00 UTC May 13 → tfiles[8]='13_12' t=0  idt=252
TARGET_TIMES = [
    {'idt_1km': 210, 'fi_111': 1, 't_111': 0},
    {'idt_1km': 241, 'fi_111': 6, 't_111': 1},
    {'idt_1km': 252, 'fi_111': 8, 't_111': 0},
]

SIM_NAMES     = ['1 km', '111 m smth-1km', '111 m']
PANEL_LABELS  = [['a', 'b', 'c'],
                 ['d', 'e', 'f'],
                 ['g', 'h', 'i']]

# ── Map domain ───────────────────────────────────────────────
latmin, latmax = 34.380, 34.530
lonmin, lonmax = -119.830, -119.570

# ── Froude computation domain ─────────────────────────────────
fr_latmin, fr_latmax = 34.38, 34.54
fr_lonmin, fr_lonmax = -119.83, -119.57

# ── Contour / level settings ─────────────────────────────────
terLvls  = np.arange(100, 1200, 300)
wspdLevs = np.arange(0, 20, 2)

# ── Tick locations ───────────────────────────────────────────
dx = np.arange(-119.8, -119.6, 0.1)
dy = np.arange(34.4,   34.5,   0.1)

# ── Lidar station locations ───────────────────────────────────
lidar_lons = [-119.6593, -119.7087, -119.7694]   # MFS-2, SBBG
lidar_lats = [ 34.4403,   34.4551,    34.4516]

# ── Constants / markers ──────────────────────────────────────
g          = 9.81
STAR_COLOR = 'white'
STAR_SIZE  = 100

OUTDIR = os.path.dirname(os.path.abspath(__file__))


# ============================================================
# ── HELPERS (mirrors existing notebook) ─────────────────────
# ============================================================
def pbl_height_thv_plus_k(thv_col, z_agl_col, thv_crit=1.5):
    ind_min = np.nanargmin(thv_col)
    val_min = thv_col[ind_min]
    above   = np.where(thv_col[ind_min:] > val_min + thv_crit)[0]
    if len(above) == 0:
        return np.nan
    ind_cross = above[0] + ind_min
    thv_lo = thv_col[ind_cross - 1]
    thv_hi = thv_col[ind_cross]
    z_lo   = z_agl_col[ind_cross - 1]
    z_hi   = z_agl_col[ind_cross]
    return ((val_min + thv_crit - thv_lo) / (thv_hi - thv_lo)) * (z_hi - z_lo) + z_lo


def froude_column(z_col, thv_col, wspd_col, zi, offset, nz):
    if zi is None or np.isnan(zi) or zi <= 0:
        zi = 1.0
    above = np.where(z_col > zi)[0]
    if len(above) == 0:
        return np.nan
    ind_zi   = max(above[0] - 1, 0)
    mean_thv = np.nanmean(thv_col[: ind_zi + 1])
    ind_free = min(ind_zi + offset, nz - 1)
    thv_free = thv_col[ind_free]
    mean_spd = np.nanmean(wspd_col[: ind_zi + 1])
    g_prime  = g * (thv_free - mean_thv) / mean_thv
    if g_prime <= 0:
        return np.nan
    return mean_spd / np.sqrt(g_prime * zi)


def compute_froude_map(ncfile, tidx, use_wrf_pblh=False):
    z   = getvar(ncfile, "z",      tidx)
    th  = getvar(ncfile, "th",     tidx)
    qv  = getvar(ncfile, "QVAPOR", tidx)
    V   = getvar(ncfile, "va",     tidx)
    ter = getvar(ncfile, "ter",    -1)
    lat2d, lon2d = latlon_coords(z)

    thv  = th * (1.0 + 0.61 * qv)
    wspd = np.abs(to_np(V))

    z_np   = to_np(z)
    thv_np = to_np(thv)
    ter_np = to_np(ter)
    lat_np = to_np(lat2d)
    lon_np = to_np(lon2d)

    nz, ny, nx = z_np.shape
    if   nz == 54:  offset = 2
    elif nz == 80:  offset = 3
    elif nz == 106: offset = 4
    else:           offset = 2

    ix_cols        = np.where((lon_np[0, :] >= fr_lonmin) &
                               (lon_np[0, :] <= fr_lonmax))[0]
    super2sub_lats = np.full(len(ix_cols), np.nan)
    lons_out       = lon_np[0, ix_cols]

    pblh2d = to_np(getvar(ncfile, "PBLH", tidx)) if use_wrf_pblh else None

    for ci, ix in enumerate(ix_cols):
        iy_rows = np.where((lat_np[:, ix] >= fr_latmin) &
                           (lat_np[:, ix] <= fr_latmax))[0]
        if len(iy_rows) < 2:
            continue

        fr_col = np.full(len(iy_rows), np.nan)
        for ri, iy in enumerate(iy_rows):
            z_vert    = z_np[:, iy, ix]
            thv_vert  = thv_np[:, iy, ix]
            wspd_vert = wspd[:, iy, ix]
            z_agl     = z_vert - ter_np[iy, ix]
            zi = pblh2d[iy, ix] if use_wrf_pblh else pbl_height_thv_plus_k(thv_vert, z_agl)
            fr_col[ri] = froude_column(z_vert, thv_vert, wspd_vert, zi, offset, nz)

        lats_col = lat_np[iy_rows, ix]
        ter_col  = ter_np[iy_rows, ix]
        if lats_col[0] < lats_col[-1]:
            lats_col = lats_col[::-1]
            fr_col   = fr_col[::-1]
            ter_col  = ter_col[::-1]

        ind_peak = np.argmax(ter_col)
        is_super = fr_col > 1.0

        for ri in range(ind_peak, len(fr_col) - 1):
            if np.isnan(fr_col[ri]) or np.isnan(fr_col[ri + 1]):
                continue
            fr_a, fr_b   = fr_col[ri],   fr_col[ri + 1]
            lat_a, lat_b = lats_col[ri], lats_col[ri + 1]
            if abs(fr_b - fr_a) < 1e-10:
                continue
            lat_cross = lat_a + (1.0 - fr_a) * (lat_b - lat_a) / (fr_b - fr_a)
            if is_super[ri] and not is_super[ri + 1]:
                if np.isnan(super2sub_lats[ci]):
                    super2sub_lats[ci] = lat_cross
                    break

    return lons_out, super2sub_lats


# ============================================================
# ── HELPER : draw one panel ──────────────────────────────────
# ============================================================
def draw_wind_panel(ax, ncfile, tidx, label, sim_name,
                    lons_fr, super2sub_lats,
                    show_ylabels=True, show_xlabels=True,
                    show_legend=False):

    timeWRF  = getvar(ncfile, 'Times',   tidx)
    U10, V10 = getvar(ncfile, 'uvmet10', tidx)
    ter      = getvar(ncfile, 'ter',     -1)
    wspd10   = np.sqrt(to_np(U10)**2 + to_np(V10)**2)
    lats2d, lons2d = latlon_coords(U10)

    ax.set_extent([lonmin, lonmax, latmin, latmax], crs=crs.PlateCarree())
    ax.set_xticks(dx, crs=crs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.set_yticks(dy, crs=crs.PlateCarree())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.tick_params(direction='out', labelsize=14, length=5, pad=2)

    if not show_ylabels:
        ax.set_yticklabels([])
    if not show_xlabels:
        ax.set_xticklabels([])

    if show_ylabels:
        ax.set_ylabel('Latitude [N]', fontsize=16)
    if show_xlabels:
        ax.set_xlabel('Longitude [W]', fontsize=16)

    states = NaturalEarthFeature(category='cultural', scale='10m',
                                 facecolor='none',
                                 name='admin_1_states_provinces')
    ax.add_feature(states, linewidth=0.5, edgecolor='black')
    ax.coastlines('10m', linewidth=2)
    ax.gridlines(color='grey', linestyle='dotted')

    wc = ax.contourf(to_np(lons2d), to_np(lats2d), wspd10,
                     levels=wspdLevs,
                     transform=crs.PlateCarree(),
                     cmap='jet', alpha=0.6, extend='max')

    ax.contour(to_np(lons2d), to_np(lats2d), to_np(ter),
               terLvls, colors='black', linewidths=0.8,
               transform=crs.PlateCarree())

    Q = ax.quiver(to_np(lons2d), to_np(lats2d),
                  to_np(U10), to_np(V10),
                  transform=crs.PlateCarree(),
                  color='k', regrid_shape=20,
                  pivot='middle', width=0.004)
    ax.quiverkey(Q, 0.78, 0.93, 10, '10 m s$^{-1}$',
                 labelpos='E', coordinates='axes',
                 fontproperties=dict(size=12))

    ax.scatter(lidar_lons, lidar_lats,
               marker='^', linewidth=3.0,
               facecolors='none', edgecolors='r', s=200,
               transform=crs.PlateCarree(), zorder=10)

    valid = ~np.isnan(super2sub_lats)
    ax.scatter(lons_fr[valid], super2sub_lats[valid],
               marker='*', color=STAR_COLOR,
               edgecolors='k', linewidths=0.4, s=STAR_SIZE,
               transform=crs.PlateCarree(), zorder=11,
               label='super→sub (Fr=1)' if show_legend else '_nolegend_')

    if show_legend:
        ax.legend(fontsize=11, loc='lower left', framealpha=0.8)

    # PDT time from WRF Times variable
    ts    = pd.to_datetime(to_np(timeWRF))
    t_pdt = (ts - pd.DateOffset(hours=7)).strftime('%H%M')

    ax.text(0.02, 0.97,
            f'({label}) {sim_name}  {t_pdt} PDT',
            transform=ax.transAxes,
            fontsize=14, va='top',
            bbox=dict(facecolor='white', edgecolor='none',
                      boxstyle='round,pad=0.2'),
            zorder=12)

    return wc


# ============================================================
# ── OPEN FILES ───────────────────────────────────────────────
# ============================================================
print('Opening 1 km file …')
ncfile_1km = Dataset(fName_1km)

nc_sm  = []
nc_111 = []
for tt in TARGET_TIMES:
    fp_sm  = baseFolder_sm  + wrffi + tfiles[tt['fi_111']] + ':00:00'
    fp_111 = baseFolder_111 + wrffi + tfiles[tt['fi_111']] + ':00:00'
    print(f'Opening {fp_sm}')
    nc_sm.append(Dataset(fp_sm))
    print(f'Opening {fp_111}')
    nc_111.append(Dataset(fp_111))

# Cartopy projection per simulation row (fixed for all time steps)
proj_1km = get_cartopy(getvar(ncfile_1km, 'uvmet10', TARGET_TIMES[0]['idt_1km']))
proj_sm  = get_cartopy(getvar(nc_sm[0],   'uvmet10', TARGET_TIMES[0]['t_111']))
proj_111 = get_cartopy(getvar(nc_111[0],  'uvmet10', TARGET_TIMES[0]['t_111']))
PROJS    = [proj_1km, proj_sm, proj_111]


# ============================================================
# ── COMPUTE FROUDE TRANSITION MAPS ───────────────────────────
# ============================================================
print('\nComputing Froude maps …')
fr_data = {}   # (row, col) → (lons_fr, super2sub_lats)

for col, tt in enumerate(TARGET_TIMES):
    print(f'  [1km]  col={col}  idt={tt["idt_1km"]}')
    fr_data[(0, col)] = compute_froude_map(ncfile_1km, tt['idt_1km'], use_wrf_pblh=False)
    print(f'  [smth] col={col}  t={tt["t_111"]}')
    fr_data[(1, col)] = compute_froude_map(nc_sm[col],  tt['t_111'], use_wrf_pblh=False)
    print(f'  [111m] col={col}  t={tt["t_111"]}')
    fr_data[(2, col)] = compute_froude_map(nc_111[col], tt['t_111'], use_wrf_pblh=False)


# ============================================================
# ── FIGURE : 3 rows × 3 cols + colorbar row ─────────────────
# ============================================================
print('\nBuilding figure …')
fig = plt.figure(figsize=(24, 20))
gs  = gridspec.GridSpec(
    4, 3,
    figure=fig,
    hspace=0.06,
    wspace=0.05,
    height_ratios=[1, 1, 1, 0.05],
)

# 3×3 subplot axes
axes = []
for row in range(3):
    row_axes = []
    for col in range(3):
        ax = fig.add_subplot(gs[row, col], projection=PROJS[row])
        row_axes.append(ax)
    axes.append(row_axes)

ax_cb = fig.add_subplot(gs[3, :])

# ncfiles and tidxs per (row, col)
ncfiles = [
    [ncfile_1km,  ncfile_1km,  ncfile_1km ],   # row 0 : 1 km
    [nc_sm[0],    nc_sm[1],    nc_sm[2]   ],   # row 1 : smth
    [nc_111[0],   nc_111[1],   nc_111[2]  ],   # row 2 : 111 m
]
tidxs = [
    [tt['idt_1km'] for tt in TARGET_TIMES],    # row 0
    [tt['t_111']   for tt in TARGET_TIMES],    # row 1
    [tt['t_111']   for tt in TARGET_TIMES],    # row 2
]

for row in range(3):
    for col in range(3):
        lons_fr, s2b = fr_data[(row, col)]
        draw_wind_panel(
            axes[row][col],
            ncfiles[row][col],
            tidxs[row][col],
            PANEL_LABELS[row][col],
            SIM_NAMES[row],
            lons_fr, s2b,
            show_ylabels = (col == 0),
            show_xlabels = (row == 2),
            show_legend  = (row == 0 and col == 0),
        )

# ── Shared colorbar ───────────────────────────────────────────
norm_cb = mcolors.BoundaryNorm(boundaries=wspdLevs, ncolors=plt.cm.jet.N)
sm_cb   = plt.cm.ScalarMappable(cmap='jet', norm=norm_cb)
sm_cb.set_array([])
cb = plt.colorbar(sm_cb, cax=ax_cb, orientation='horizontal', extend='max')
cb.set_label('10-m wind speed (m s$^{-1}$)', fontsize=17)
cb.ax.tick_params(labelsize=15)
cb.set_ticks(wspdLevs)

# ── Save ─────────────────────────────────────────────────────
fbase = 'SWEX2022_IOP10_2Dwinds_threeTimes_thv1p5K'
for ext, dpi in [('jpg', 150), ('pdf', 300)]:
    fig.savefig(os.path.join(OUTDIR, f'{fbase}.{ext}'), bbox_inches='tight', dpi=dpi)
    print(f'Saved {fbase}.{ext}')

plt.close()

# ── Close files ───────────────────────────────────────────────
ncfile_1km.close()
for nc in nc_sm + nc_111:
    nc.close()

print('Done.')
