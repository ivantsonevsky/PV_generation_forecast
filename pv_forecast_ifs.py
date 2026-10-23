#!/usr/bin/env python3
"""
PV output forecast from the ECMWF IFS (Control) sunshine duration forecast.

Data
  IFS Control from MARS (via the ECMWF Web API, ~/.ecmwfapirc):
    sund (param 189) - sunshine duration, accumulated from the start of the run [s]
    2t   (param 167) - 2 m temperature [K], used for the module temperature
    lcc, mcc, hcc (params 186/187/188) - low, medium and high cloud cover [0-1]
    tcc  (param 164) - total cloud cover [0-1], used to set the layer overlap
  Steps are hourly up to +90 h (Control forecast) has hourly output there) and 3-hourly
  after that, so a 3-day forecast is fully hourly.

Method (per forecast interval, evaluated every 5 minutes inside it)
  1. Sun position and Ineichen-Perez clear-sky irradiance (GHI, DNI, DHI).
  2. Sunshine fraction f = IFS sunshine seconds / seconds in which the sun
     *could* shine. WMO sunshine means DNI > 120 W/m2, so only the minutes with
     clear-sky DNI above 120 W/m2 count as possible sunshine.
  3. All-sky irradiance: beam = f * DNI_cs,
     diffuse = f * DHI_cs + (1 - f) * T_cloud * GHI_cs
     T_cloud is the light that gets through the clouds hiding the sun. Each
     layer passes a fraction of the clear-sky GHI (low ~25 %, medium ~45 %,
     high/cirrus ~75 %); with random overlap of the IFS low/medium/high cloud
     cover, T_cloud is the mean transmission given that at least one layer
     covers the sun:
       T_cloud = (through - clear) / (1 - clear)
     where, under random overlap, through = prod(1 - c_i + c_i t_i) and
     clear = prod(1 - c_i). IFS total cloud cover (tcc) uses maximum-random
     overlap, so the layers are usually more stacked than random. The overlap is
     therefore a blend of random and maximum overlap, weighted so that the clear
     fraction equals 1 - tcc (stacked layers: more clear sky, thicker cloud,
     lower T_cloud). tcc only sets the overlap; the layer covers set the
     transmission. Use --no-tcc for plain random overlap.
     With almost no cloud it uses the cover-weighted mean of t_i, and with none
     at all it falls back to k_overcast (0.3).
  4. Plane-of-array irradiance for each array (isotropic sky, ground albedo,
     incidence-angle losses), NOCT cell temperature, PVWatts-style DC power,
     system losses, inverter efficiency and optional inverter AC limit. The
     arrays share one inverter: when the total is clipped, every array is
     scaled down by the same factor.

Arrays
  Any number of arrays (up to 7), each with its own capacity [kWp], tilt from
  horizontal [deg] and azimuth clockwise from north [deg] (90=E, 180=S, 270=W).
  Set them in ARRAYS below, or on the command line with --arrays N and one
  value per array for --capacity/--tilt/--azimuth/--name (a single value is
  used for all arrays). The plot shows each array and the total.

Examples
  python3 pv_forecast_ifs.py                                   # ARRAYS below
  python3 pv_forecast_ifs.py --lat 51.42 --lon -0.95 --tz Europe/London \
      --capacity 4.2 --tilt 35 --azimuth 180
  python3 pv_forecast_ifs.py --arrays 3 --capacity 4 3 2 --tilt 30 30 15 \
      --azimuth 180 90 270 --name South East West --inverter-kw 6
  # same, one array per option: kWp,tilt,azimuth[,name]
  python3 pv_forecast_ifs.py --array 4,30,180,South --array 3,30,90,East --array 2,15,270,West
"""
import argparse
import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# ---- defaults: change these to your site and system ----
LAT, LON, ELEVATION_M = 51.42, -0.95, 73.0
TIMEZONE = "Europe/London"
# one entry per array: name, capacity [kWp], tilt [deg], azimuth [deg from north]
ARRAYS = [
    {"name": "SW", "kwp": 3.5, "tilt": 30.0, "azimuth": 230.0},
    # {"name": "East", "kwp": 3.0, "tilt": 30.0, "azimuth": 90.0},
    # {"name": "West", "kwp": 2.0, "tilt": 15.0, "azimuth": 270.0},
]
DEFAULT_TILT, DEFAULT_AZIMUTH = 30.0, 180.0
MAX_ARRAYS = 7              # one colour each in the plot

MARS_LATENCY_H = 7          # Control FC 00/12 UTC runs are in MARS ~6-7 h after base time
HOURLY_MAX_STEP = 90        # Control FC has hourly steps up to +90 h, 3-hourly up to +144 h
SUBSTEP_MIN = 5
SUNSHINE_DNI = 120.0        # WMO sunshine threshold [W/m2]
CACHE_DIR = Path(__file__).resolve().parent / "ifs_cache"


# ----------------------------------------------------------------- IFS data
def candidate_runs(now_utc):
    """Most recent Control FC runs that should already be in MARS, newest first."""
    t = now_utc - dt.timedelta(hours=MARS_LATENCY_H)
    base = t.replace(hour=12 if t.hour >= 12 else 0, minute=0, second=0, microsecond=0)
    return [base, base - dt.timedelta(hours=12)]


def mars_steps(max_step):
    steps = list(range(0, min(max_step, HOURLY_MAX_STEP) + 1))
    steps += list(range(HOURLY_MAX_STEP + 3, max_step + 3, 3))
    return steps


def retrieve_ifs(base, lat, lon, max_step):
    """Download sund, 2t, tcc and the cloud layers for a small box around the
    site (cached). Returns the valid times and a dict of point time series."""
    from ecmwfapi import ECMWFService

    CACHE_DIR.mkdir(exist_ok=True)
    target = CACHE_DIR / f"ifs_{base:%Y%m%d%H}_{lat:.2f}_{lon:.2f}_{max_step}_tcc.nc"
    if not target.exists():
        request = {
            "class": "od", "stream": "oper", "type": "fc", "expver": "1",
            "levtype": "sfc", "param": "189/167/186/187/188/164",
            "date": f"{base:%Y-%m-%d}", "time": f"{base:%H}",
            "step": "/".join(map(str, mars_steps(max_step))),
            "area": f"{lat + 0.2}/{lon - 0.2}/{lat - 0.2}/{lon + 0.2}",
            "grid": "0.1/0.1", "format": "netcdf",
        }
        tmp = target.with_suffix(".part")
        ECMWFService("mars", verbose=False).execute(request, str(tmp))
        tmp.rename(target)
    with xr.open_dataset(target, engine="scipy") as ds:
        p = ds.isel(latitude=int(abs(ds.latitude - lat).argmin()),
                    longitude=int(abs(ds.longitude - lon).argmin())).load()
    times = pd.DatetimeIndex(p.time.values).tz_localize("UTC")
    ifs = {v: p[v].values.astype(float) for v in ("sund", "lcc", "mcc", "hcc", "tcc")}
    ifs["t2m"] = p.t2m.values.astype(float) - 273.15
    return times, ifs


# --------------------------------------------------------- solar geometry
def solar_position(times, lat, lon):
    """Zenith and azimuth (deg, azimuth clockwise from north) and the
    extraterrestrial normal irradiance, NOAA approximations."""
    doy = times.dayofyear.values
    hour = times.hour.values + times.minute.values / 60 + times.second.values / 3600
    g = 2 * np.pi / 365 * (doy - 1 + (hour - 12) / 24)
    eqtime = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g)
                       - 0.014615 * np.cos(2 * g) - 0.040849 * np.sin(2 * g))
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g)
            - 0.006758 * np.cos(2 * g) + 0.000907 * np.sin(2 * g)
            - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    ha = np.radians((hour * 60 + eqtime + 4 * lon) / 4 - 180)
    phi = np.radians(lat)
    cosz = np.sin(phi) * np.sin(decl) + np.cos(phi) * np.cos(decl) * np.cos(ha)
    zenith = np.degrees(np.arccos(np.clip(cosz, -1, 1)))
    azimuth = np.degrees(np.arctan2(np.sin(ha), np.cos(ha) * np.sin(phi)
                                    - np.tan(decl) * np.cos(phi))) + 180
    e0 = (1.000110 + 0.034221 * np.cos(g) + 0.001280 * np.sin(g)
          + 0.000719 * np.cos(2 * g) + 0.000077 * np.sin(2 * g))
    return zenith, azimuth, 1361.0 * e0


def clearsky_ineichen(zenith, dni_extra, elevation, linke):
    """Ineichen-Perez clear-sky GHI, DNI, DHI [W/m2]."""
    cosz = np.cos(np.radians(zenith))
    up = zenith < 90
    am = np.full_like(zenith, np.nan)
    am[up] = 1 / (cosz[up] + 0.50572 * (96.07995 - zenith[up]) ** -1.6364)
    am_abs = am * np.exp(-elevation / 8434.5)
    fh1, fh2 = np.exp(-elevation / 8000), np.exp(-elevation / 1250)
    cg1, cg2 = 5.09e-5 * elevation + 0.868, 3.92e-5 * elevation + 0.0387
    ghi = cg1 * dni_extra * cosz * np.exp(-cg2 * am_abs * (fh1 + fh2 * (linke - 1)))
    dni = dni_extra * (0.664 + 0.163 / fh1) * np.exp(-0.09 * am_abs * (linke - 1))
    dni2 = ghi * (1 - (0.1 - 0.2 * np.exp(-linke)) / (0.1 + 0.882 / fh1)) / np.where(up, cosz, 1)
    dni = np.minimum(dni, dni2)
    ghi, dni = np.where(up, ghi, 0.0), np.where(up, dni, 0.0)
    dhi = np.maximum(ghi - dni * cosz, 0.0)
    return ghi, dni, dhi


def plane_of_array(dni, dhi, ghi, zenith, sun_az, tilt, azimuth, albedo):
    z, b = np.radians(zenith), np.radians(tilt)
    cos_aoi = np.cos(z) * np.cos(b) + np.sin(z) * np.sin(b) * np.cos(np.radians(sun_az - azimuth))
    cos_aoi = np.clip(cos_aoi, 0, 1)
    iam = np.clip(1 - 0.05 * (1 / np.maximum(cos_aoi, 1e-3) - 1), 0, 1)  # ASHRAE, b0=0.05
    return (dni * cos_aoi * iam
            + dhi * (1 + np.cos(b)) / 2
            + ghi * albedo * (1 - np.cos(b)) / 2)


# ------------------------------------------------------------- PV model
CLOUD_TRANSMISSION = {"lcc": 0.25, "mcc": 0.45, "hcc": 0.75}


def cloud_transmission(cover, k_overcast, tcc=None):
    """Mean fraction of the clear-sky GHI that reaches the ground while the sun
    is hidden. The layer covers set the transmission; tcc (if given) sets how
    much the layers overlap (random overlap if tcc is None)."""
    layers = list(CLOUD_TRANSMISSION)
    c = np.stack([np.clip(cover[l], 0, 1) for l in layers]).astype(float)
    t = np.array([CLOUD_TRANSMISSION[l] for l in layers], dtype=float)[:, None] * np.ones_like(c)

    through = np.prod(1 - c + c * t, axis=0)      # random overlap
    clear = np.prod(1 - c, axis=0)

    if tcc is not None:
        # maximum overlap: layers nested, largest cover first
        order = np.argsort(-c, axis=0)
        cs, ts = np.take_along_axis(c, order, 0), np.take_along_axis(t, order, 0)
        nxt = np.vstack([cs[1:], np.zeros((1, c.shape[1]))])
        clear_m = 1 - cs[0]
        through_m = clear_m + ((cs - nxt) * np.cumprod(ts, axis=0)).sum(axis=0)
        # weight w (0 = random, 1 = maximum) so that the clear fraction is 1 - tcc
        span = clear_m - clear
        w = np.where(span > 1e-6,
                     np.clip((1 - np.clip(tcc, 0, 1) - clear) / np.maximum(span, 1e-6), 0, 1),
                     0.0)
        through = through + w * (through_m - through)
        clear = clear + w * span

    with np.errstate(invalid="ignore", divide="ignore"):
        trans = (through - clear) / (1 - clear)
        tot = c.sum(axis=0)
        thin = np.where(tot > 1e-9, (c * t).sum(axis=0) / tot, k_overcast)
    return np.where(1 - clear > 1e-3, trans, thin)


def pv_forecast(times, ifs, arrays, a):
    """5-minute AC power [kW] for the IFS sunshine forecast and for clear sky."""
    t0 = times[:-1]
    dt_s = np.diff(times.asi8) / 1e9
    sun_s = np.clip(np.diff(ifs["sund"]), 0, None)

    sub, idx = [], []
    for i, (start, length) in enumerate(zip(t0, dt_s)):
        n = int(round(length / 60 / SUBSTEP_MIN))
        sub.append(start + pd.to_timedelta((np.arange(n) + 0.5) * SUBSTEP_MIN, unit="min"))
        idx.append(np.full(n, i))
    sub, idx = sub[0].append(sub[1:]), np.concatenate(idx)

    zen, saz, dni_extra = solar_position(sub, a.lat, a.lon)
    ghi_cs, dni_cs, dhi_cs = clearsky_ineichen(zen, dni_extra, a.elevation, a.linke)

    # sunshine fraction relative to the possible sunshine in each interval
    possible = np.bincount(idx, weights=(dni_cs > SUNSHINE_DNI) * SUBSTEP_MIN * 60.0,
                           minlength=len(t0))
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(possible > 0, np.clip(sun_s / possible, 0, 1), np.nan)
    frac = pd.Series(frac).interpolate(limit_direction="both").fillna(0).values
    f = frac[idx]

    # cloud cover and temperature are instantaneous: interpolate in time
    at_sub = lambda v: np.interp(sub.asi8, times.asi8, v)  # noqa: E731
    cover = {layer: at_sub(ifs[layer]) for layer in CLOUD_TRANSMISSION}
    tcc = None if a.no_tcc else np.clip(at_sub(ifs["tcc"]), 0, 1)
    t_cloud = cloud_transmission(cover, a.k_overcast, tcc)
    tair = at_sub(ifs["t2m"])

    dni = f * dni_cs
    dhi = f * dhi_cs + (1 - f) * t_cloud * ghi_cs
    ghi = dni * np.cos(np.radians(zen)) + dhi

    def ac_power(dni, dhi, ghi):
        """AC power of each array and of the system [kW]."""
        per = []
        for arr in arrays:
            poa = plane_of_array(dni, dhi, ghi, zen, saz, arr["tilt"], arr["azimuth"], a.albedo)
            tcell = tair + poa / 800 * (a.noct - 20)
            dc = np.maximum(arr["kwp"] * poa / 1000 * (1 + a.gamma * (tcell - 25)), 0)
            per.append(dc * (1 - a.losses) * a.inverter_eff)
        per = np.array(per)
        total = per.sum(axis=0)
        if a.inverter_kw:  # shared inverter: clip the total, scale arrays alike
            scale = np.minimum(1, a.inverter_kw / np.maximum(total, 1e-9))
            per, total = per * scale, total * scale
        return per, total

    per, total = ac_power(dni, dhi, ghi)
    _, total_cs = ac_power(dni_cs, dhi_cs, ghi_cs)
    return pd.DataFrame({
        "pv_kw": total,
        **{f"arr{i}_kw": per[i] for i in range(len(arrays))},
        "clearsky_kw": total_cs,
        "sunshine_frac": f,
        "ghi_wm2": ghi,
        "t2m_c": tair,
        **{f"{layer}_pct": 100 * cover[layer] for layer in CLOUD_TRANSMISSION},
        "cloud_trans": t_cloud,
    }, index=sub)


# ------------------------------------------------------------------ plot
def plot(hourly, daily, base, a, arrays, out):
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    ink, ink2, grid, surface = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    blue, gray = "#2a78d6", "#8a8984"   # total, clear sky
    colors = ["#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
    multi = len(arrays) > 1
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": grid, "axes.labelcolor": ink2,
                         "xtick.color": ink2, "ytick.color": ink2, "text.color": ink})

    fig = plt.figure(figsize=(11, 9), facecolor=surface)
    outer = fig.add_gridspec(2, 1, height_ratios=[4.2, 2], hspace=0.42)
    top = outer[0].subgridspec(2, 1, height_ratios=[3, 1.1], hspace=0.38)
    ax1 = fig.add_subplot(top[0])
    axc = fig.add_subplot(top[1], sharex=ax1)
    ax2 = fig.add_subplot(outer[1])
    for ax in (ax1, axc, ax2):
        ax.set_facecolor(surface)
        ax.grid(axis="y", color=grid, linewidth=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=0)

    x = hourly.index.tz_localize(None) + pd.Timedelta(minutes=30)
    ax1.plot(x, hourly.clearsky_kw, color=gray, lw=1.5, ls=(0, (4, 3)),
             label="Clear-sky potential (total)" if multi else "Clear-sky potential")
    ax1.fill_between(x, hourly.pv_kw, color=blue, alpha=0.12, lw=0)
    ax1.plot(x, hourly.pv_kw, color=blue, lw=2.4,
             label="Total" if multi else "Forecast (IFS sunshine duration)")
    if multi:
        for i, arr in enumerate(arrays):
            ax1.plot(x, hourly[f"arr{i}_kw"], color=colors[i], lw=1.6,
                     label=f"{arr['name']}  {describe(arr)}")
    for d in daily.index[1:]:
        ax1.axvline(pd.Timestamp(d), color=grid, lw=1)
    ax1.set_ylabel("AC power [kW]")
    ax1.set_ylim(0, max(hourly.clearsky_kw.max(), 0.1) * 1.12)
    ax1.set_xlim(pd.Timestamp(daily.index[0]), pd.Timestamp(daily.index[-1]) + pd.Timedelta(days=1))
    ax1.tick_params(which="both", labelbottom=False)

    # cloud cover strips, sky order: high on top, low at the bottom
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("cloud", [surface, "#b9b8b2", "#3d3c39"])
    t_edges = mdates.date2num(hourly.index.tz_localize(None).append(
        hourly.index[-1:].tz_localize(None) + pd.Timedelta(hours=1)))
    cover = np.vstack([hourly.lcc_pct, hourly.mcc_pct, hourly.hcc_pct])
    mesh = axc.pcolormesh(t_edges, [0, 1, 2, 3], cover, cmap=cmap, vmin=0, vmax=100,
                          shading="flat", edgecolor=surface, linewidth=0.3)
    for d in daily.index[1:]:
        axc.axvline(pd.Timestamp(d), color=ink2, lw=0.8, alpha=0.4)
    axc.set_yticks([0.5, 1.5, 2.5])
    axc.set_yticklabels(["Low", "Medium", "High"])
    axc.set_ylim(0, 3)
    axc.grid(False)
    axc.set_title("IFS cloud cover", loc="left", fontweight="bold", fontsize=10)
    cax = axc.inset_axes([0.82, 1.1, 0.18, 0.09])
    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal", ticks=[0, 50, 100])
    cb.ax.set_xticklabels(["0%", "50%", "100%"])
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=0, labelsize=8)
    cb.ax.xaxis.set_ticks_position("top")
    axc.xaxis.set_major_locator(mdates.HourLocator(byhour=12))
    axc.xaxis.set_major_formatter(mdates.DateFormatter("12:00\n%a %d %b"))
    axc.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 6, 18]))
    axc.xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
    axc.tick_params(which="both", length=0, labelcolor=ink2)
    ncol = 3 if multi else 2
    rows = int(np.ceil((len(arrays) + 2 if multi else 2) / ncol))
    ax1.legend(loc="lower left", frameon=False, ncol=ncol, bbox_to_anchor=(0, 1.0),
               borderaxespad=0.2, handlelength=2.2)
    ax1.set_title("Hourly PV output", loc="left", fontweight="bold", pad=8 + 19 * rows)

    labels = [f"{pd.Timestamp(d):%a %d %b}\n{r.sunshine_h:.1f} h sunshine"
              for d, r in daily.iterrows()]
    pos = np.arange(len(daily))
    # one bar per array, then the total (with the clear-sky potential outlined)
    nbars = len(arrays) + 1 if multi else 1
    slot = 0.84 / nbars if multi else 0.55
    off = lambda j: (j - (nbars - 1) / 2) * slot  # noqa: E731
    if multi:
        for i, arr in enumerate(arrays):
            bars = ax2.bar(pos + off(i), daily[f"arr{i}_kwh"], width=slot * 0.88,
                           color=colors[i], label=arr["name"])
            for b, v in zip(bars, daily[f"arr{i}_kwh"]):
                ax2.text(b.get_x() + b.get_width() / 2, v, f"{v:.1f}",
                         ha="center", va="bottom", fontsize=8, color=ink2)
    xt = pos + off(nbars - 1)
    ax2.bar(xt, daily.clearsky_kwh, width=slot * 0.88, color="none", edgecolor=gray,
            lw=1.2, ls=(0, (4, 3)), label="Clear-sky potential")
    bars = ax2.bar(xt, daily.energy_kwh, width=slot * 0.88, color=blue,
                   label="Total" if multi else "Forecast")
    for b, v, top in zip(bars, daily.energy_kwh, np.maximum(daily.energy_kwh, daily.clearsky_kwh)):
        ax2.text(b.get_x() + b.get_width() / 2, top, f"{v:.1f} kWh" if not multi else f"{v:.1f}",
                 ha="center", va="bottom", fontsize=10 if not multi else 9,
                 color=ink, fontweight="bold")
    ax2.set_xticks(pos)
    ax2.set_xticklabels(labels)
    ax2.set_ylabel("Energy [kWh/day]")
    ax2.set_ylim(0, max(daily.clearsky_kwh.max(), 0.1) * 1.18)
    ax2.set_title("Daily PV energy", loc="left", fontweight="bold")
    ax2.legend(loc="lower right", frameon=False, ncol=nbars + 1, bbox_to_anchor=(1, 1.0),
               borderaxespad=0.2)

    arr_txt = (f"{len(arrays)} arrays, {sum(r['kwp'] for r in arrays):g} kWp" if multi
               else describe(arrays[0]))
    loc = f"{abs(a.lat):.2f}°{'N' if a.lat >= 0 else 'S'} {abs(a.lon):.2f}°{'E' if a.lon >= 0 else 'W'}"
    fig.suptitle(f"PV forecast  {loc}  |  {arr_txt}",
                 x=0.06, ha="left", fontsize=13, fontweight="bold")
    fig.text(0.06, 0.935, f"ECMWF IFS Control {base:%Y-%m-%d %H} UTC run, sunshine duration and low/medium/high "
             f"cloud cover (hourly steps); "
             f"times in {a.tz}", fontsize=9, color=ink2)
    fig.subplots_adjust(left=0.075, right=0.98, top=0.87 - 0.045 * (rows - 1), bottom=0.08)
    fig.savefig(out, dpi=130, facecolor=surface)
    return fig


# ------------------------------------------------------------------ main
def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lat", type=float, default=LAT)
    p.add_argument("--lon", type=float, default=LON)
    p.add_argument("--elevation", type=float, default=ELEVATION_M, help="site height [m]")
    p.add_argument("--tz", default=TIMEZONE, help="time zone for days and plot (IANA name)")
    p.add_argument("--arrays", type=int, metavar="N",
                   help=f"number of PV arrays (1-{MAX_ARRAYS}); default: ARRAYS in the script")
    p.add_argument("--capacity", type=float, nargs="+", metavar="KWP",
                   help="array size [kWp], one value per array")
    p.add_argument("--tilt", type=float, nargs="+", metavar="DEG",
                   help=f"panel tilt from horizontal [deg], one per array or one for all (default {DEFAULT_TILT:g})")
    p.add_argument("--azimuth", type=float, nargs="+", metavar="DEG",
                   help="panel facing, clockwise from north [deg]: 90=E, 180=S, 270=W; "
                        f"one per array or one for all (default {DEFAULT_AZIMUTH:g})")
    p.add_argument("--name", nargs="+", help="array names, one per array (default: Array 1, 2, ...)")
    p.add_argument("--array", action="append", metavar="KWP,TILT,AZ[,NAME]",
                   help="one array; repeat for each (instead of --arrays/--capacity/...)")
    p.add_argument("--inverter-kw", type=float, default=None, help="inverter AC limit [kW]")
    p.add_argument("--days", type=int, default=3, help="forecast days, 1-5 (hourly data covers 3)")
    p.add_argument("--run", help="IFS run 'YYYY-MM-DD HH' (UTC); default: latest available")
    p.add_argument("--losses", type=float, default=0.10, help="DC system losses (wiring, soiling, mismatch)")
    p.add_argument("--inverter-eff", type=float, default=0.96)
    p.add_argument("--gamma", type=float, default=-0.0037, help="power temperature coefficient [1/K]")
    p.add_argument("--noct", type=float, default=45.0, help="nominal operating cell temperature [C]")
    p.add_argument("--albedo", type=float, default=0.2)
    p.add_argument("--linke", type=float, default=3.0, help="Linke turbidity for clear sky (2 clean - 5 hazy)")
    p.add_argument("--k-overcast", type=float, default=0.3,
                   help="overcast GHI as fraction of clear-sky GHI when IFS has no cloud in any layer")
    p.add_argument("--no-tcc", action="store_true",
                   help="ignore IFS total cloud cover and assume random overlap of the layers")
    p.add_argument("--out", default="pv_forecast.png", help="plot file")
    p.add_argument("--csv", default="pv_forecast_hourly.csv", help="hourly output table")
    p.add_argument("--show", action="store_true", help="also open the plot window")
    return p.parse_args()


def build_arrays(a):
    """Array list from --array, from --arrays/--capacity/..., or from ARRAYS."""
    if a.array:
        arrays = []
        for i, spec in enumerate(a.array, 1):
            v = spec.split(",")
            if len(v) not in (3, 4):
                raise SystemExit(f"--array {spec}: expected KWP,TILT,AZ[,NAME]")
            arrays.append({"name": v[3] if len(v) == 4 else f"Array {i}",
                           "kwp": float(v[0]), "tilt": float(v[1]), "azimuth": float(v[2])})
    elif a.arrays or a.capacity:
        n = a.arrays or len(a.capacity)

        def per_array(values, default, option):
            if values is None:
                if default is None:
                    raise SystemExit(f"--{option} is needed, one value per array ({n})")
                values = [default]
            if len(values) == 1:
                values = values * n
            if len(values) != n:
                raise SystemExit(f"--{option}: give 1 or {n} values (one per array), got {len(values)}")
            return values

        names = per_array(a.name, None, "name") if a.name else [f"Array {i}" for i in range(1, n + 1)]
        arrays = [{"name": nm, "kwp": k, "tilt": t, "azimuth": z} for nm, k, t, z in zip(
            names, per_array(a.capacity, None, "capacity"),
            per_array(a.tilt, DEFAULT_TILT, "tilt"), per_array(a.azimuth, DEFAULT_AZIMUTH, "azimuth"))]
    else:
        arrays = [dict(arr) for arr in ARRAYS]
    if not 1 <= len(arrays) <= MAX_ARRAYS:
        raise SystemExit(f"Give between 1 and {MAX_ARRAYS} arrays, not {len(arrays)}.")
    return arrays


def describe(arr):
    return f"{arr['kwp']:g} kWp {arr['tilt']:g}°/{arr['azimuth']:g}°"


def main():
    a = parse_args()
    if not a.show:
        import matplotlib
        matplotlib.use("Agg")
    arrays = build_arrays(a)
    tz = a.tz

    runs = ([dt.datetime.strptime(a.run, "%Y-%m-%d %H").replace(tzinfo=dt.timezone.utc)]
            if a.run else candidate_runs(dt.datetime.now(dt.timezone.utc)))
    for base in runs:
        # first local day whose sunrise the run covers
        base_local = pd.Timestamp(base).tz_convert(tz)
        first = base_local.normalize()
        if base_local > first + pd.Timedelta(hours=6):
            first += pd.Timedelta(days=1)
        end = first + pd.Timedelta(days=a.days)
        max_step = int(np.ceil((end - pd.Timestamp(base)) / pd.Timedelta(hours=1)))
        max_step = min(144, max_step + (-max_step) % (3 if max_step > HOURLY_MAX_STEP else 1))
        try:
            print(f"IFS Control {base:%Y-%m-%d %H} UTC, steps 0-{max_step} h ...")
            times, ifs = retrieve_ifs(base, a.lat, a.lon, max_step)
            break
        except Exception as e:  # run not yet in MARS -> try the previous one
            print(f"  not available ({e})")
    else:
        raise SystemExit("No IFS run could be retrieved.")

    df = pv_forecast(times, ifs, arrays, a)
    sund = ifs["sund"]
    df.index = df.index.tz_convert(tz)
    df = df[(df.index >= first) & (df.index < end)]

    hourly = df.drop(columns="sunshine_frac").resample("1h").mean()
    # IFS sunshine minutes per hour, straight from the accumulated field
    edges = hourly.index.append(hourly.index[-1:] + pd.Timedelta(hours=1)).asi8
    hourly["sunshine_min"] = np.diff(np.interp(edges, times.asi8, sund)) / 60
    names = {f"arr{i}_kw": f"{arr['name']}_kw" for i, arr in enumerate(arrays)}
    hourly.rename(columns=names).round(3).to_csv(a.csv, index_label=f"time ({a.tz})")

    step_h = SUBSTEP_MIN / 60
    daily = pd.DataFrame({
        "energy_kwh": df.pv_kw.resample("1D").sum() * step_h,
        "clearsky_kwh": df.clearsky_kw.resample("1D").sum() * step_h,
        "peak_kw": df.pv_kw.resample("1D").max(),
        **{f"arr{i}_kwh": df[f"arr{i}_kw"].resample("1D").sum() * step_h
           for i in range(len(arrays))},
    })
    # mean cloud cover of each layer while the sun is up
    day = df[df.clearsky_kw > 0]
    for layer in CLOUD_TRANSMISSION:
        daily[f"{layer}_pct"] = day[f"{layer}_pct"].resample("1D").mean()
    # IFS sunshine hours per local day, straight from the accumulated field
    edges = [d.value for d in daily.index] + [(daily.index[-1] + pd.Timedelta(days=1)).value]
    daily["sunshine_h"] = np.diff(np.interp(edges, times.asi8, sund)) / 3600
    daily.index = daily.index.date

    print("\nPV arrays:")
    for arr in arrays:
        print(f"  {arr['name']}: {describe(arr)}")
    table = daily.rename(columns={f"arr{i}_kwh": f"{arr['name']}_kwh" for i, arr in enumerate(arrays)})
    print(table.round(2).to_string())
    print(f"\nTotal {daily.energy_kwh.sum():.1f} kWh over {len(daily)} days"
          + "".join(f"; {arr['name']} {daily[f'arr{i}_kwh'].sum():.1f}"
                    for i, arr in enumerate(arrays) if len(arrays) > 1))

    plot(hourly, daily, base, a, arrays, a.out)
    print(f"Plot: {Path(a.out).resolve()}\nHourly table: {Path(a.csv).resolve()}")
    if a.show:
        import matplotlib.pyplot as plt
        plt.show()


if __name__ == "__main__":
    main()
