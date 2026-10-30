#!/usr/bin/env python
"""
Station time series near the SBBG and MFS2 lidars during SWEX IOP10, to relate
surface wind reversals (rotor / subcritical flow) to temperature and dew point.

For each lidar site the nearest SCE, RAWS and ISFS stations are selected and a
figure is made with T, Td, wind direction, wind speed (+ gusts) and a station map.
ISFS gusts are the max 3-s running-mean speed per 5 min from the 20 Hz sonic
winds (ISFS/v10/high_rate). RAWS gusts (hourly, as reported) come from the older
Mesowest download data_full/obs/Stations_SB_20220401_20220515_PDT_lvl1.mat; the
SCE/RAWS Synoptic netCDF files have no gust variable. Each
observation is also classified into a flow regime (used for the transitions CSV):
  downslope : northerly flow,  v < -V_THRESHOLD   (jet reaches the surface)
  reversed  : southerly flow,  v > +V_THRESHOLD   (rotor / onshore return flow)
  weak      : |v| <= V_THRESHOLD
Regime transitions that persist for at least --persist minutes are listed in a
CSV with the change in T, Td and RH (mean over --window minutes after minus
before the transition).

Loaders are shared with plot_stations_IOP10.py (same folder).

Usage:
  python timeseries_rotor_stations_IOP10.py
  python timeseries_rotor_stations_IOP10.py --radius 4 --nmax 8
  python timeseries_rotor_stations_IOP10.py --stations SE504 SE741 351SE s04
"""

import argparse
import os

import numpy as np
import pandas as pd
import xarray as xr
from scipy.io import loadmat
import cartopy.crs as ccrs
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from plot_stations_IOP10 import (DATA_ROOT, DOMAINS, LIDARS, MS_TO_KNOTS, NETWORK_STYLE,
                                 PDT_OFFSET, draw_basemap, load_isfs, load_synoptic,
                                 load_topography)

plt.rcParams.update({'font.size': 14})   # publication (supplement) sizing

V_THRESHOLD = 0.5   # m/s, |v| below this is "weak"


def distance_km(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.deg2rad, (lon1, lat1, lon2, lat2))
    a = (np.sin((lat2 - lat1) / 2) ** 2 +
         np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2)
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def flow_regime(df):
    """+1 downslope (northerly), -1 reversed (southerly), 0 weak, NaN missing."""
    v_ms = df['v'] / MS_TO_KNOTS
    reg = pd.Series(np.where(v_ms < -V_THRESHOLD, 1, np.where(v_ms > V_THRESHOLD, -1, 0)),
                    index=df.index, dtype=float)
    return reg.where(v_ms.notna())


def find_transitions(st, persist, window):
    """Persistent downslope <-> reversed transitions (weak periods are bridged)."""
    df = st['data']
    reg = flow_regime(df)
    reg_f = reg.replace(0, np.nan).ffill()
    prev = reg_f.shift()
    change = reg_f.notna() & prev.notna() & (reg_f != prev)

    rows = []
    for t in reg_f.index[change.values]:
        new = reg_f[t]
        after = reg[(reg.index >= t) & (reg.index < t + persist)].replace(0, np.nan).dropna()
        if after.empty or (after != new).any():
            continue
        before_df = df[(df.index >= t - window) & (df.index < t)]
        after_df = df[(df.index >= t) & (df.index < t + window)]
        rows.append({
            'station': st['stid'], 'network': st['network'], 'elev_m': st['elev_m'],
            'time_utc': t, 'time_pdt': t - PDT_OFFSET,
            'transition': 'N->S (reversal)' if new == -1 else 'S->N (downslope)',
            'dT': after_df['T'].mean() - before_df['T'].mean(),
            'dTd': after_df['Td'].mean() - before_df['Td'].mean(),
            'dRH': after_df['RH'].mean() - before_df['RH'].mean(),
            'spd_before_ms': before_df['spd'].mean(),
            'spd_after_ms': after_df['spd'].mean(),
        })
    return rows, reg


def load_isfs_gusts(folder, nums, t0, t1, gust_s=3, interval='5min'):
    """Gusts (max of the gust_s running-mean speed per interval) from the hourly
    20 Hz ISFS files isfs_swex_qc_hr_YYYYMMDD_HH.nc, U/V_10m_sN(time, sample)."""
    pieces = {n: [] for n in nums}
    nsamp = 20
    for h in pd.date_range(t0.floor('H'), t1.floor('H'), freq='H'):
        f = os.path.join(folder, f'isfs_swex_qc_hr_{h:%Y%m%d_%H}.nc')
        if not os.path.exists(f):
            print('missing', f)
            continue
        with xr.open_dataset(f) as ds:
            nsamp = ds.sizes['sample']
            t = (ds['time'].values[:, None] +
                 (np.arange(nsamp) * 1e9 / nsamp).astype('timedelta64[ns]')).ravel()
            for n in nums:
                if f'U_10m_s{n}' not in ds:
                    continue
                spd = np.hypot(ds[f'U_10m_s{n}'].values, ds[f'V_10m_s{n}'].values).ravel()
                pieces[n].append(pd.Series(spd, index=t))
    gusts = {}
    for n, p in pieces.items():
        if not p:
            continue
        s = pd.concat(p).sort_index()
        s = s[~s.index.duplicated()]
        g = s.rolling(f'{gust_s}s', min_periods=int(0.75 * gust_s * nsamp)).mean()
        g = g.resample(interval).max()
        # ISFS 5-min means are stamped at the interval centre; do the same here
        g.index = g.index + pd.Timedelta(interval) / 2
        gusts[n] = g
    print(f'ISFS gusts computed for towers: {sorted(gusts, key=int)}')
    return gusts


def load_mat_gusts(path):
    """Hourly gusts per station ID from the Mesowest lvl1 .mat (48 stations x hours,
    MATLAB datenum in PDT, NaN-filled) -> {stid: Series indexed in UTC}."""
    m = loadmat(path)
    names = [str(np.squeeze(s)) for s in m['stations'].ravel()]
    tnum, gust = m['obs_tnum_PDT_lvl1'], m['obs_wspd_gust_lvl1']
    out = {}
    for i, stid in enumerate(names):
        ok = tnum[i] > 0
        t = pd.to_datetime(tnum[i, ok] - 719529, unit='D').round('s') + PDT_OFFSET
        g = pd.Series(gust[i, ok], index=t)
        g = g[~g.index.duplicated()].sort_index()
        if g.notna().any():
            out[stid] = g
    print(f'gusts from {os.path.basename(path)}: {len(out)} stations')
    return out


def select_stations(stations, site, radius_km, nmax, explicit):
    name, lon, lat = site
    # shallow copies, so each lidar keeps its own dist_km (data frames are shared)
    stations = [dict(st, dist_km=distance_km(lon, lat, st['lon'], st['lat'])) for st in stations]
    if explicit:
        sel = [st for st in stations if st['stid'] in explicit]
    else:
        sel = [st for st in stations if st['dist_km'] <= radius_km]
    sel = sorted(sel, key=lambda s: s['dist_km'])[:nmax]
    print(f'{name}: ' + ', '.join(f"{s['stid']} ({s['dist_km']:.1f} km, {s['elev_m']:.0f} m)"
                                  for s in sel))
    return sel


def plot_station_map(ax, site, sel, colors, topo, show_site=True):
    """Terrain map around the lidar with the selected stations in their line colors.
    show_site=False for station groups without a lidar (site is then the group centre)."""
    name, lon0, lat0 = site
    dmax = max([st['dist_km'] for st in sel] + [1.0]) + 0.8   # km, margin around stations
    dlat = dmax / 111.0
    dlon = dlat / np.cos(np.deg2rad(lat0))
    extent = [lon0 - dlon, lon0 + dlon, lat0 - dlat, lat0 + dlat]
    draw_basemap(ax, extent, topo)
    ax.tick_params(labelsize=11)

    if show_site:
        ax.scatter(lon0, lat0, s=110, marker='s', color='black', zorder=10)
        t = ax.annotate(name, xy=(lon0, lat0), xytext=(0, -10), textcoords='offset points',
                        ha='center', va='top', fontsize=13, fontweight='bold', zorder=11)
        t.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),
                            path_effects.Normal()])

    for i, st in enumerate(sel):
        style = NETWORK_STYLE[st['network']]
        ax.scatter(st['lon'], st['lat'], s=100, marker=style['marker'], color=colors[i % 10],
                   edgecolor='black', linewidth=0.8, zorder=10)
        # put the label on the side facing away from the lidar so it does not cover it
        dx = (st['lon'] - lon0) * np.cos(np.deg2rad(lat0))
        dy = st['lat'] - lat0
        ang = np.arctan2(dy, dx) if (dx or dy) else np.pi / 4
        ox, oy = 7 * np.cos(ang), 7 * np.sin(ang)
        t = ax.annotate(st['stid'], xy=(st['lon'], st['lat']), xytext=(ox, oy),
                        textcoords='offset points', fontsize=12, fontweight='bold',
                        ha='left' if ox >= 0 else 'right', va='bottom' if oy >= 0 else 'top',
                        color=colors[i % 10], zorder=11)
        t.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),
                            path_effects.Normal()])

    handles = [Line2D([], [], ls='none', marker=s['marker'], markerfacecolor='white',
                      markeredgecolor='black', label=s['label']) for s in NETWORK_STYLE.values()]
    ax.legend(handles=handles, fontsize=11, loc='lower left', framealpha=0.9)


def plot_site(site, sel, t0, t1, highlight, transitions, topo, outdir, tag='IOP10',
              exts=('png', 'pdf'), show_site=True):
    name = site[0]
    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    fig = plt.figure(figsize=(19, 13))
    gs = fig.add_gridspec(4, 2, width_ratios=[2.3, 1], wspace=0.12, hspace=0.22)
    ax_t = fig.add_subplot(gs[0, 0])
    axs = [ax_t] + [fig.add_subplot(gs[k, 0], sharex=ax_t) for k in range(1, 4)]
    ax_t, ax_td, ax_dir, ax_spd = axs
    for ax in axs[:-1]:
        ax.tick_params(labelbottom=False)
    # right column: station legend (row 0, next to T), map (rows 1-3)
    ax_map = fig.add_subplot(gs[1:4, 1], projection=ccrs.PlateCarree())
    plot_station_map(ax_map, site, sel, colors, topo, show_site)
    for ax, letter in zip(axs + [ax_map], 'abcde'):
        ax.text(0.0, 1.02, f'({letter})', transform=ax.transAxes, ha='left', va='bottom',
                fontsize=15, fontweight='bold')

    has_gust = False
    for i, st in enumerate(sel):
        df = st['data']
        df = df[(df.index >= t0) & (df.index <= t1)]
        tp = df.index - PDT_OFFSET
        c = colors[i % 10]
        label = f"{st['stid']} ({st['network']}, {st['elev_m']:.0f} m, {st['dist_km']:.1f} km)"
        ax_t.plot(tp, df['T'], color=c, lw=1.3, label=label)
        ax_td.plot(tp, df['Td'], color=c, lw=1.3)
        ax_dir.scatter(tp, df['dir'], color=c, s=8)
        ax_spd.plot(tp, df['spd'], color=c, lw=1.3)
        gust = st.get('gust', df['gust'] if 'gust' in df else None)
        if gust is not None:
            gust = gust[(gust.index >= t0) & (gust.index <= t1)]
            ax_spd.plot(gust.index - PDT_OFFSET, gust, color=c, lw=1.0, ls='--')
            has_gust = True

    if has_gust:
        ax_spd.legend(handles=[Line2D([], [], color='black', lw=1.3, label='mean'),
                               Line2D([], [], color='black', lw=1.0, ls='--',
                                      label='gust')],
                      fontsize=11, loc='upper right', framealpha=0.9)

    for ax in axs:
        ax.axvspan(highlight[0] - PDT_OFFSET, highlight[1] - PDT_OFFSET,
                   color='khaki', alpha=0.3, zorder=0)
        ax.grid(alpha=0.3)

    ax_t.set_ylabel('T (°C)')
    ax_td.set_ylabel('T$_d$ (°C)')
    ax_dir.set_ylabel('Wind dir.')
    ax_dir.set_ylim(0, 360)
    ax_dir.set_yticks([0, 90, 180, 270, 360])
    ax_dir.set_yticklabels(['N', 'E', 'S', 'W', 'N'])
    ax_spd.set_ylabel('Wind speed (m s$^{-1}$)')
    # explicit PDT ticks every 2 h with string labels (mdates locators misplace ticks
    # with the pandas datetime converter in this env)
    t0_pdt, t1_pdt = t0 - PDT_OFFSET, t1 - PDT_OFFSET
    ticks = pd.date_range(t0_pdt.ceil('2H'), t1_pdt, freq='2H')
    ax_spd.set_xlim(t0_pdt, t1_pdt)
    ax_spd.set_xticks(ticks)
    ax_spd.set_xticklabels([t.strftime('%H:%M\n%d %b') for t in ticks])
    ax_spd.set_xlabel('Time (PDT)')

    ax_t.legend(fontsize=12, loc='center left', bbox_to_anchor=(1.02, 0.5),
                title='Station (network, elev., dist. to ' + ('lidar)' if show_site else 'centre)'),
                title_fontsize=12)

    for ext in exts:
        fname = os.path.join(outdir, f'timeseries_rotor_{name}_{tag}.{ext}')
        fig.savefig(fname, bbox_inches='tight', dpi=300)
        print('saved', fname)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-root', default=DATA_ROOT)
    p.add_argument('--start', default='2022-05-12 16:00', help='PDT (IOP10 start)')
    p.add_argument('--end', default='2022-05-13 07:00', help='PDT (IOP10 end)')
    p.add_argument('--highlight', nargs=2, default=['2022-05-12 22:00', '2022-05-13 01:00'],
                   help='PDT window shaded in the figures')
    p.add_argument('--radius', type=float, default=3.0, help='km around each lidar')
    p.add_argument('--nmax', type=int, default=6, help='max stations per lidar')
    p.add_argument('--stations', nargs='*', default=None, help='explicit station IDs')
    p.add_argument('--persist', type=int, default=20, help='min. minutes a new regime must last')
    p.add_argument('--window', type=int, default=30, help='minutes averaged before/after')
    p.add_argument('--gust-seconds', type=int, default=3, help='gust averaging time (s)')
    p.add_argument('--no-gust', action='store_true', help='skip reading the ISFS high-rate files')
    p.add_argument('--gust-mat', default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'data_full', 'obs',
        'Stations_SB_20220401_20220515_PDT_lvl1.mat'), help='Mesowest .mat with RAWS gusts')
    p.add_argument('--outdir', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures'))
    args = p.parse_args()

    t0 = pd.Timestamp(args.start) + PDT_OFFSET      # -> UTC
    t1 = pd.Timestamp(args.end) + PDT_OFFSET
    highlight = [pd.Timestamp(h) + PDT_OFFSET for h in args.highlight]
    persist = pd.Timedelta(minutes=args.persist)
    window = pd.Timedelta(minutes=args.window)
    os.makedirs(args.outdir, exist_ok=True)

    root = args.data_root
    stations = (load_synoptic(os.path.join(root, 'SCE'), 'SCE') +
                load_synoptic(os.path.join(root, 'RAWS'), 'RAWS') +
                load_isfs(os.path.join(root, 'ISFS', 'v10_single_station')))
    topo = load_topography(os.path.join(root, 'geographic_files', 'sb_wrf_topo_50m.nc'),
                           DOMAINS['lee'])

    sels = [select_stations(stations, site, args.radius, args.nmax, args.stations)
            for site in LIDARS]
    if not args.no_gust:
        nums = sorted({st['num'] for sel in sels for st in sel if st['network'] == 'ISFS'})
        gusts = load_isfs_gusts(os.path.join(root, 'ISFS', 'v10', 'high_rate'), nums,
                                t0, t1, gust_s=args.gust_seconds)
        if os.path.exists(args.gust_mat):
            gusts.update(load_mat_gusts(args.gust_mat))
        else:
            print('no RAWS gusts, file not found:', args.gust_mat)
        for sel in sels:
            for st in sel:
                key = st['num'] if st['network'] == 'ISFS' else st['stid']
                if st['network'] in ('ISFS', 'RAWS') and key in gusts:
                    st['gust'] = gusts[key]

    all_rows = []
    for site, sel in zip(LIDARS, sels):
        rows = []
        for st in sel:
            tr, _ = find_transitions(st, persist, window)
            rows += [dict(r, lidar=site[0], dist_km=st['dist_km'])
                     for r in tr if t0 <= r['time_utc'] <= t1]
        plot_site(site, sel, t0, t1, highlight, rows, topo, args.outdir)
        all_rows += rows

    out = pd.DataFrame(all_rows)
    if not out.empty:
        out = out.sort_values(['lidar', 'time_utc'])
        fname = os.path.join(args.outdir, 'rotor_transitions_IOP10.csv')
        out.to_csv(fname, index=False, float_format='%.2f')
        print('saved', fname)
        rev = out[out['transition'].str.startswith('N->S')]
        if not rev.empty:
            print(f'\nN->S reversals: n={len(rev)}, median dT={rev["dT"].median():.2f} °C, '
                  f'median dTd={rev["dTd"].median():.2f} °C')
    else:
        print('no persistent regime transitions found')


if __name__ == '__main__':
    main()
