#!/usr/bin/env python3
"""
stellar_geolocate.py
────────────────────────────────────────────────────────────────────────────────
Stellar Geolocation Tool  v3.0
Determines a photographer's geographic coordinates from a night-sky image
using plate-solve data and corrected spherical astronomy.

Plate-solve pipeline (tried in order):
  1. nova.astrometry.net  — cloud API; uploads image, returns full WCS
  2. local solve-field    — astrometry.net offline binary (same algorithm)
  3. manual entry         — user types RA/Dec per star from any solver

Geolocation pipeline:
  - Multi-star least-squares optimization (scipy)
  - Correct LHA formula: LHA = GMST + lon - RA
  - Bennett atmospheric refraction correction
  - Multiple optimizer starting points to avoid local minima
  - Single-star meridian-transit fallback when scipy unavailable

Dependencies:
  pip install requests   needed for nova API and local-solver WCS fetch
  pip install scipy      needed for multi-star optimization (strongly recommended)
  pip install astropy    optional; improves WCS pixel-to-RA/Dec accuracy

Usage:
  python stellar_geolocate.py
────────────────────────────────────────────────────────────────────────────────
"""

import math
import os
import sys
import json
import time
import platform
import subprocess
import tempfile
import shutil
from datetime import datetime, timezone

# ── Optional dependencies ──────────────────────────────────────────────────────

try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False

try:
    from scipy.optimize import minimize
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

try:
    from astropy.wcs import WCS as AstropyWCS
    from astropy.io.fits import getheader as fits_getheader
    import io as _io
    ASTROPY_AVAILABLE = True
except ImportError:
    ASTROPY_AVAILABLE = False


# ══════════════════════════════════════════════════════════════════════════════
#  CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

NOVA_API_BASE     = "http://nova.astrometry.net/api"
NOVA_UPLOAD_URL   = "http://nova.astrometry.net/api/upload"
NOVA_WCS_URL      = "http://nova.astrometry.net/wcs-file/{job_id}"
POLL_INTERVAL_SEC = 8
POLL_TIMEOUT_SEC  = 360
REQUEST_TIMEOUT   = 20
GRID_COLS         = 5
GRID_ROWS         = 4


# ══════════════════════════════════════════════════════════════════════════════
#  DISPLAY HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def banner():
    print()
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║          Stellar Geolocation Tool  v3.0                         ║")
    print("║   Geolocate a photographer from a night-sky image               ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    print()

    warnings = []
    if not REQUESTS_AVAILABLE:
        warnings.append(("requests", "nova API and local-solve WCS disabled"))
    if not SCIPY_AVAILABLE:
        warnings.append(("scipy",    "multi-star optimization unavailable"))
    if not ASTROPY_AVAILABLE:
        warnings.append(("astropy",  "using gnomonic WCS approximation"))

    if warnings:
        width = max(len(f"{p}: {r}") for p, r in warnings) + 2
        print(f"  ┌─ OPTIONAL DEPENDENCIES {'─' * (width - 20)}┐")
        for pkg, reason in warnings:
            line = f"{pkg}: {reason}"
            print(f"  │  ✗  {line:<{width}}│")
        fix = "pip install " + " ".join(p for p, _ in warnings)
        print(f"  │  Fix: {fix:<{width - 2}}│")
        print(f"  └{'─' * (width + 4)}┘")
        print()


def section(title):
    print(f"\n── {title} {'─' * max(0, 62 - len(title))}")


def progress(msg, end=False):
    if end:
        print(f"\r    {msg:<60}")
    else:
        print(f"\r    {msg:<60}", end="", flush=True)


# ══════════════════════════════════════════════════════════════════════════════
#  INPUT HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def prompt_float(label, min_val=None, max_val=None, default=None):
    while True:
        suffix = f"  (default: {default})" if default is not None else ""
        raw = input(f"    {label}{suffix}: ").strip()
        if raw == "" and default is not None:
            return float(default)
        try:
            val = float(raw)
        except ValueError:
            print("      x Enter a numeric value.")
            continue
        if min_val is not None and val < min_val:
            print(f"      x Minimum is {min_val}.")
            continue
        if max_val is not None and val > max_val:
            print(f"      x Maximum is {max_val}.")
            continue
        return val


def prompt_int(label, min_val=None, max_val=None, default=None):
    while True:
        suffix = f"  (default: {default})" if default is not None else ""
        raw = input(f"    {label}{suffix}: ").strip()
        if raw == "" and default is not None:
            return int(default)
        try:
            val = int(raw)
        except ValueError:
            print("      x Enter a whole number.")
            continue
        if min_val is not None and val < min_val:
            print(f"      x Minimum is {min_val}.")
            continue
        if max_val is not None and val > max_val:
            print(f"      x Maximum is {max_val}.")
            continue
        return val


def prompt_str(label, allow_empty=False):
    while True:
        raw = input(f"    {label}: ").strip()
        if raw or allow_empty:
            return raw
        print("      x This field cannot be empty.")


def prompt_choice(label, options):
    for i, (display, _) in enumerate(options, 1):
        print(f"    {i}) {display}")
    while True:
        raw = input(f"    {label}: ").strip()
        try:
            idx = int(raw)
            if 1 <= idx <= len(options):
                return options[idx - 1][1]
            print(f"      x Enter 1-{len(options)}.")
        except ValueError:
            print("      x Enter a number.")


def prompt_yes_no(label, default=True):
    hint = "[Y/n]" if default else "[y/N]"
    raw = input(f"    {label} {hint}: ").strip().lower()
    if raw == "":
        return default
    return raw in ("y", "yes")


# ══════════════════════════════════════════════════════════════════════════════
#  MAC TIMESTAMP EXTRACTION
# ══════════════════════════════════════════════════════════════════════════════

def get_mac_timestamps(filepath):
    """
    Return Modified, Accessed, and (if available) Created timestamps as UTC datetimes.

    Platform notes:
      Windows : st_ctime = creation time
      macOS   : st_birthtime = actual birth time
      Linux   : st_ctime = inode change time (NOT creation; ext4 does not expose it
                without debugfs)
    """
    stat     = os.stat(filepath)
    modified = datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc)
    accessed = datetime.fromtimestamp(stat.st_atime, tz=timezone.utc)

    sys_name = platform.system()
    if sys_name == "Windows":
        created = datetime.fromtimestamp(stat.st_ctime, tz=timezone.utc)
        created_label = "Created (Windows)"
    elif sys_name == "Darwin":
        try:
            created = datetime.fromtimestamp(stat.st_birthtime, tz=timezone.utc)
            created_label = "Created (macOS birthtime)"
        except AttributeError:
            created, created_label = None, None
    else:
        created = datetime.fromtimestamp(stat.st_ctime, tz=timezone.utc)
        created_label = "Inode Change (Linux -- NOT creation time)"

    return {
        "modified": (modified, "Modified"),
        "accessed": (accessed, "Accessed"),
        "created":  (created, created_label) if created else None,
    }


def display_timestamps(timestamps, filepath):
    print(f"\n    File : {filepath}")
    for key in ("modified", "accessed", "created"):
        entry = timestamps.get(key)
        if entry:
            dt, label = entry
            print(f"    {label:<35} {dt.strftime('%Y-%m-%d  %H:%M:%S  UTC')}")


# ══════════════════════════════════════════════════════════════════════════════
#  DATE / TIME COLLECTION
# ══════════════════════════════════════════════════════════════════════════════

def prompt_datetime_manual():
    print()
    print("    Enter observation date and time in UTC.")
    print("    Convert from local time before entering if necessary.")
    print()
    while True:
        year   = prompt_int("Year",                  min_val=1900, max_val=2100)
        month  = prompt_int("Month          (1-12)", min_val=1,    max_val=12)
        day    = prompt_int("Day            (1-31)", min_val=1,    max_val=31)
        hour   = prompt_int("Hour  UTC      (0-23)", min_val=0,    max_val=23)
        minute = prompt_int("Minute         (0-59)", min_val=0,    max_val=59)
        second = prompt_int("Second         (0-59)", min_val=0,    max_val=59)
        try:
            return datetime(year, month, day, hour, minute, second, tzinfo=timezone.utc)
        except ValueError as e:
            print(f"\n      x Invalid date: {e}. Please re-enter.\n")


def collect_observation_datetime(image_path):
    section("Step 2 -- Observation Date & Time")
    timestamps = None
    if image_path:
        try:
            timestamps = get_mac_timestamps(image_path)
            print()
            print("  MAC timestamps detected:")
            display_timestamps(timestamps, image_path)
        except Exception as exc:
            print(f"\n    [WARNING] Could not read timestamps: {exc}")

    if timestamps:
        print()
        print("  Which timestamp best represents the capture time?")
        print("  Modified is most commonly the capture time for camera transfers.")
        print()
        options = []
        for key in ("modified", "accessed", "created"):
            entry = timestamps.get(key)
            if entry:
                dt, label = entry
                options.append((
                    f"{label:<35} {dt.strftime('%Y-%m-%d  %H:%M:%S  UTC')}",
                    dt
                ))
        options.append(("Enter date and time manually", None))
        chosen = prompt_choice("Choice", options)
        if chosen is not None:
            return chosen

    return prompt_datetime_manual()


# ══════════════════════════════════════════════════════════════════════════════
#  WCS CONVERSION
# ══════════════════════════════════════════════════════════════════════════════

def gnomonic_to_radec(xi_deg, eta_deg, ra0_deg, dec0_deg):
    """
    Full gnomonic (TAN) inverse projection: tangent-plane offsets to RA/Dec.
    xi_deg  : east  offset in degrees on tangent plane
    eta_deg : north offset in degrees on tangent plane
    Returns (ra_deg, dec_deg).
    """
    xi   = math.radians(xi_deg)
    eta  = math.radians(eta_deg)
    ra0  = math.radians(ra0_deg)
    dec0 = math.radians(dec0_deg)

    rho = math.sqrt(xi ** 2 + eta ** 2)
    if rho < 1e-14:
        return ra0_deg, dec0_deg

    c     = math.atan(rho)
    sin_c = math.sin(c)
    cos_c = math.cos(c)

    dec = math.asin(
        cos_c * math.sin(dec0) + eta * sin_c * math.cos(dec0) / rho
    )
    ra = ra0 + math.atan2(
        xi * sin_c,
        rho * math.cos(dec0) * cos_c - eta * math.sin(dec0) * sin_c
    )

    return math.degrees(ra) % 360.0, math.degrees(dec)


def pixel_to_radec_gnomonic(px, py, cx, cy, pixscale_arcsec, orientation_deg,
                             ra0_deg, dec0_deg, parity=1.0):
    """
    Convert image pixel coordinates to RA/Dec using astrometry.net calibration
    and the full gnomonic projection. No FITS or astropy required.

    px, py           : pixel coordinates (top-left origin, +x right, +y down)
    cx, cy           : image center pixel
    pixscale_arcsec  : plate scale in arcseconds per pixel
    orientation_deg  : position angle of image +Y axis, measured N toward E
    parity           : +1 standard astronomical (east = left when N up),
                       -1 mirrored. astrometry.net parity=1 means standard.
    """
    s  = pixscale_arcsec / 3600.0
    pa = math.radians(orientation_deg)

    # Pixel offset from center; flip y so +eta is up (north when PA=0)
    xi_pix  =  (px - cx)
    eta_pix = -(py - cy)

    # Parity: standard astronomical has east to the left (negative x in image)
    xi_sky_pix  = -parity * xi_pix
    eta_sky_pix =  eta_pix

    # Scale to degrees
    xi_deg  = xi_sky_pix  * s
    eta_deg = eta_sky_pix * s

    # Rotate by position angle
    d_east  =  xi_deg * math.cos(pa) + eta_deg * math.sin(pa)
    d_north = -xi_deg * math.sin(pa) + eta_deg * math.cos(pa)

    return gnomonic_to_radec(d_east, d_north, ra0_deg, dec0_deg)


def pixel_to_radec_astropy(px, py, wcs_obj):
    """Convert pixel to RA/Dec using astropy WCS. Returns (ra_deg, dec_deg)."""
    sky = wcs_obj.all_pix2world([[px, py]], 0)[0]
    return float(sky[0]) % 360.0, float(sky[1])


def build_astropy_wcs(fits_bytes):
    """Parse FITS bytes with astropy. Returns a WCS object or None."""
    if not ASTROPY_AVAILABLE:
        return None
    try:
        buf = _io.BytesIO(fits_bytes)
        hdr = fits_getheader(buf, ignore_missing_simple=True)
        return AstropyWCS(hdr)
    except Exception:
        return None


def sample_stars_grid(image_params, calib, wcs_obj=None):
    """
    Sample a GRID_COLS x GRID_ROWS grid of sky positions across the image.
    Returns list of star dicts with ra_deg, dec_deg, alt_apparent, alt_corrected.
    """
    W  = image_params["image_width"]
    H  = image_params["image_height"]
    cx = W / 2.0
    cy = H / 2.0

    stars = []
    for row in range(GRID_ROWS):
        for col in range(GRID_COLS):
            px = W * (0.1 + 0.8 * col / (GRID_COLS - 1))
            py = H * (0.1 + 0.8 * row / (GRID_ROWS - 1))

            if wcs_obj is not None:
                ra_deg, dec_deg = pixel_to_radec_astropy(px, py, wcs_obj)
            else:
                ra_deg, dec_deg = pixel_to_radec_gnomonic(
                    px, py, cx, cy,
                    calib["pixscale"],
                    calib["orientation"],
                    calib["ra"],
                    calib["dec"],
                    calib.get("parity", 1.0)
                )

            alt_apparent  = _altitude_from_pixel(py, image_params)
            refraction    = _refraction_correction(alt_apparent)
            alt_corrected = alt_apparent - refraction

            stars.append({
                "ra_deg":        ra_deg,
                "dec_deg":       dec_deg,
                "pixel_y":       py,
                "alt_apparent":  alt_apparent,
                "alt_corrected": alt_corrected,
                "source":        "wcs_grid",
            })

    above = [s for s in stars if s["alt_apparent"] >= 0]
    return above if above else stars


# ══════════════════════════════════════════════════════════════════════════════
#  NOVA API
# ══════════════════════════════════════════════════════════════════════════════

def check_nova_connectivity():
    """Returns (True, None) if reachable, (False, reason) if not."""
    if not REQUESTS_AVAILABLE:
        return False, "requests library not installed"
    try:
        resp = requests.get(NOVA_API_BASE + "/", timeout=8)
        resp.raise_for_status()
        return True, None
    except requests.exceptions.ConnectionError:
        return False, "connection refused or DNS failure"
    except requests.exceptions.Timeout:
        return False, "connection timed out"
    except requests.exceptions.HTTPError as e:
        return False, f"HTTP {e.response.status_code}"
    except Exception as e:
        return False, str(e)


def nova_login(api_key):
    resp = requests.post(
        NOVA_API_BASE + "/login",
        data={"request-json": json.dumps({"apikey": api_key})},
        timeout=REQUEST_TIMEOUT
    )
    resp.raise_for_status()
    result = resp.json()
    if result.get("status") != "success":
        raise RuntimeError(f"Login failed: {result.get('error', result)}")
    return result["session"]


def nova_upload(session, image_path):
    with open(image_path, "rb") as fh:
        resp = requests.post(
            NOVA_UPLOAD_URL,
            data={"request-json": json.dumps({
                "session":              session,
                "allow_commercial_use": "n",
                "allow_modifications":  "n",
                "publicly_visible":     "n",
            })},
            files={"file": fh},
            timeout=120,
        )
    resp.raise_for_status()
    result = resp.json()
    if result.get("status") != "success":
        raise RuntimeError(f"Upload failed: {result.get('error', result)}")
    return result["subid"]


def nova_poll_submission(subid):
    deadline = time.time() + POLL_TIMEOUT_SEC
    attempt  = 0
    while time.time() < deadline:
        attempt += 1
        try:
            resp   = requests.get(f"{NOVA_API_BASE}/submissions/{subid}", timeout=REQUEST_TIMEOUT)
            result = resp.json()
            jobs   = result.get("jobs", [])
            if jobs and jobs[0] is not None:
                return jobs[0]
        except Exception:
            pass
        elapsed = int(time.time() - (deadline - POLL_TIMEOUT_SEC))
        progress(f"Waiting for job assignment... {elapsed}s (attempt {attempt})")
        time.sleep(POLL_INTERVAL_SEC)
    raise TimeoutError("Timed out waiting for job assignment.")


def nova_poll_job(job_id):
    deadline = time.time() + POLL_TIMEOUT_SEC
    attempt  = 0
    while time.time() < deadline:
        attempt += 1
        try:
            resp   = requests.get(f"{NOVA_API_BASE}/jobs/{job_id}", timeout=REQUEST_TIMEOUT)
            status = resp.json().get("status")
            if status == "success":
                progress("Solve completed successfully.", end=True)
                return True
            if status == "failure":
                raise RuntimeError(
                    "Plate solve failed -- image may not contain enough stars "
                    "or the FOV may be unusually large or small."
                )
        except RuntimeError:
            raise
        except Exception:
            pass
        elapsed = int(time.time() - (deadline - POLL_TIMEOUT_SEC))
        progress(f"Solving... {elapsed}s (job {job_id}, attempt {attempt})")
        time.sleep(POLL_INTERVAL_SEC)
    raise TimeoutError("Timed out waiting for plate solve to complete.")


def nova_get_calibration(job_id):
    resp = requests.get(f"{NOVA_API_BASE}/jobs/{job_id}/calibration", timeout=REQUEST_TIMEOUT)
    resp.raise_for_status()
    cal = resp.json()
    required = ("ra", "dec", "radius", "pixscale", "orientation", "parity")
    missing  = [k for k in required if k not in cal]
    if missing:
        raise RuntimeError(f"Calibration missing keys: {missing}")
    return cal


def nova_get_wcs_fits(job_id):
    try:
        url  = NOVA_WCS_URL.format(job_id=job_id)
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        if len(resp.content) > 100:
            return resp.content
    except Exception:
        pass
    return None


def run_nova_solve(image_path):
    """
    Full nova.astrometry.net pipeline.
    Returns (calib_dict, wcs_bytes_or_None).
    """
    print()
    print("  A free API key is required.")
    print("  Register at https://nova.astrometry.net  then find your key under My Profile.")
    print()
    api_key = prompt_str("  API key")

    print()
    progress("Logging in...")
    session = nova_login(api_key)
    progress("Logged in. Uploading image...")
    subid = nova_upload(session, image_path)
    progress(f"Uploaded (sub {subid}). Waiting for job assignment...")
    job_id = nova_poll_submission(subid)
    progress(f"Job {job_id} assigned. Solving...")
    nova_poll_job(job_id)

    progress("Fetching calibration...")
    calib = nova_get_calibration(job_id)

    progress("Downloading WCS FITS...")
    wcs_bytes = nova_get_wcs_fits(job_id)
    if wcs_bytes:
        progress("WCS file downloaded.", end=True)
    else:
        progress("WCS unavailable; gnomonic approximation will be used.", end=True)

    return calib, wcs_bytes


# ══════════════════════════════════════════════════════════════════════════════
#  LOCAL SOLVE-FIELD
# ══════════════════════════════════════════════════════════════════════════════

def local_solver_path():
    """Return path to solve-field binary, or None."""
    found = shutil.which("solve-field")
    if found:
        return found
    for candidate in [
        "/usr/local/astrometry/bin/solve-field",
        "/opt/astrometry/bin/solve-field",
        "/usr/local/bin/solve-field",
        os.path.expanduser("~/astrometry/bin/solve-field"),
    ]:
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None


def parse_local_wcs_file(wcs_path):
    """
    Parse a .wcs FITS file from solve-field. Requires astropy.
    Returns a calib dict (matching nova format) with extra _wcs_obj key, or None.
    """
    if not ASTROPY_AVAILABLE:
        return None
    try:
        hdr  = fits_getheader(wcs_path, ignore_missing_simple=True)
        wcs  = AstropyWCS(hdr)
        ra0  = float(hdr.get("CRVAL1", 0))
        dec0 = float(hdr.get("CRVAL2", 0))

        if "CD1_1" in hdr and "CD1_2" in hdr:
            cd11            = float(hdr["CD1_1"])
            cd12            = float(hdr["CD1_2"])
            pixscale_deg    = math.sqrt(cd11 ** 2 + cd12 ** 2)
            orientation_deg = math.degrees(math.atan2(cd12, -cd11))
            parity          = -1.0 if (cd11 * float(hdr.get("CD2_2", 1)) < 0) else 1.0
        elif "CDELT1" in hdr:
            pixscale_deg    = abs(float(hdr["CDELT1"]))
            orientation_deg = float(hdr.get("CROTA2", 0.0))
            parity          = 1.0
        else:
            return None

        return {
            "ra":          ra0,
            "dec":         dec0,
            "pixscale":    pixscale_deg * 3600.0,
            "orientation": orientation_deg,
            "parity":      parity,
            "radius":      0.0,
            "_wcs_obj":    wcs,
        }
    except Exception as e:
        print(f"\n    [WARNING] Could not parse WCS: {e}")
        return None


def run_local_solver(image_path, solver_bin):
    """
    Run solve-field. Returns (calib_dict, wcs_obj_or_None) or raises.
    """
    print()
    print(f"  Solver binary : {solver_bin}")
    print()

    workdir  = tempfile.mkdtemp(prefix="stellar_geo_")
    basename = os.path.splitext(os.path.basename(image_path))[0]
    out_base = os.path.join(workdir, basename)
    wcs_path = out_base + ".wcs"

    cmd = [
        solver_bin,
        "--no-plots",
        "--overwrite",
        "--new-fits", "none",
        "--wcs",      wcs_path,
        "--out",      out_base,
        image_path,
    ]

    progress("Running solve-field...")
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        shutil.rmtree(workdir, ignore_errors=True)
        raise TimeoutError("solve-field timed out after 5 minutes.")
    except Exception as e:
        shutil.rmtree(workdir, ignore_errors=True)
        raise RuntimeError(f"solve-field could not run: {e}")

    if proc.returncode != 0 or not os.path.exists(wcs_path):
        diag = (proc.stderr or "")[-400:].strip() or "(no stderr)"
        shutil.rmtree(workdir, ignore_errors=True)
        raise RuntimeError(
            f"solve-field exited {proc.returncode}.\nDiagnostic: {diag}"
        )

    progress("Parsing WCS output...")
    calib = parse_local_wcs_file(wcs_path)
    shutil.rmtree(workdir, ignore_errors=True)

    if calib is None:
        raise RuntimeError(
            "solve-field succeeded but WCS could not be parsed. "
            "Install astropy:  pip install astropy"
        )

    wcs_obj = calib.pop("_wcs_obj", None)
    progress("Local solve complete.", end=True)
    return calib, wcs_obj


# ══════════════════════════════════════════════════════════════════════════════
#  PLATE SOLVE ORCHESTRATOR
# ══════════════════════════════════════════════════════════════════════════════

class SolveMethod:
    NOVA   = "nova.astrometry.net API"
    LOCAL  = "local solve-field"
    MANUAL = "manual entry"


def _apply_calib(calib, wcs_obj, image_params, method_label):
    """Shared post-solve: display calibration, update FOV, sample star grid."""
    print()
    print(f"  Plate center : RA {calib['ra']:.4f}  Dec {calib['dec']:+.4f}")
    print(f"  Plate scale  : {calib['pixscale']:.3f} arcsec/pixel")
    print(f"  Orientation  : {calib['orientation']:.2f} degrees")
    wcs_src = "astropy FITS WCS" if wcs_obj is not None else "gnomonic approximation"
    print(f"  WCS source   : {wcs_src}")

    computed_fov = calib["pixscale"] * image_params["image_height"] / 3600.0
    print(f"  Computed FOV : {computed_fov:.2f} degrees (vertical)")
    if abs(computed_fov - image_params.get("fov_vertical", 0)) > 0.5:
        image_params["fov_vertical"]      = computed_fov
        image_params["degrees_per_pixel"] = computed_fov / image_params["image_height"]
        print(f"  Updated FOV in image params.")

    stars = sample_stars_grid(image_params, calib, wcs_obj)
    print(f"  Sampled {len(stars)} sky positions from WCS grid ({GRID_COLS}x{GRID_ROWS}).")
    return stars


def run_plate_solve(image_path, image_params):
    """
    Attempt plate solving in priority order:
      1. nova.astrometry.net API
      2. local solve-field binary
      3. manual star entry

    Returns (stars_list, method_string).
    """
    section("Step 4 -- Plate Solve")
    attempts = []

    # ── 1. Nova ──────────────────────────────────────────────────────────────
    print()
    reachable, reason = check_nova_connectivity()
    if not reachable:
        print(f"  nova.astrometry.net : UNREACHABLE ({reason})")
    else:
        print(f"  nova.astrometry.net : reachable")
        if prompt_yes_no("  Use nova API to auto-plate-solve the image?"):
            try:
                calib, wcs_bytes = run_nova_solve(image_path)
                wcs_obj = build_astropy_wcs(wcs_bytes) if wcs_bytes else None
                stars   = _apply_calib(calib, wcs_obj, image_params, SolveMethod.NOVA)
                return stars, SolveMethod.NOVA
            except KeyboardInterrupt:
                print("\n    Cancelled.")
            except Exception as exc:
                attempts.append((SolveMethod.NOVA, str(exc)))
                print(f"\n    nova failed: {exc}")
                print("    Trying next method.")

    # ── 2. Local solver ───────────────────────────────────────────────────────
    solver_bin = local_solver_path()
    if solver_bin is None:
        print("\n  local solve-field : NOT FOUND")
        print("    Install: sudo apt install astrometry.net   (Debian/Ubuntu)")
        print("             brew install astrometry.net        (macOS)")
    else:
        print(f"\n  local solve-field : {solver_bin}")
        if prompt_yes_no("  Use local solve-field to plate-solve the image?"):
            try:
                calib, wcs_obj = run_local_solver(image_path, solver_bin)
                stars = _apply_calib(calib, wcs_obj, image_params, SolveMethod.LOCAL)
                return stars, SolveMethod.LOCAL
            except KeyboardInterrupt:
                print("\n    Cancelled.")
            except Exception as exc:
                attempts.append((SolveMethod.LOCAL, str(exc)))
                print(f"\n    Local solve failed: {exc}")
                print("    Falling back to manual entry.")

    # ── 3. Manual fallback ────────────────────────────────────────────────────
    if attempts:
        print()
        print("  ┌─ Automated solve attempts ─────────────────────────────────┐")
        for method, err in attempts:
            short = err[:36]
            print(f"  │  {method:<25}  FAILED: {short:<36}│")
        print("  └────────────────────────────────────────────────────────────┘")

    print("\n  Proceeding with manual star entry.")
    stars = collect_stars_manual(image_params)
    return stars, SolveMethod.MANUAL


# ══════════════════════════════════════════════════════════════════════════════
#  IMAGE GEOMETRY COLLECTION
# ══════════════════════════════════════════════════════════════════════════════

def collect_image_parameters():
    """
    Collect image geometry. horizon_pixel is always required for altitude
    measurements regardless of solve method.
    """
    section("Step 3 -- Image Geometry")
    print()
    print("  Y=0 is the TOP of the image.")
    print("  Horizon pixel: Y-coordinate of the visible horizon in the image.")
    print("  FOV may be updated automatically if a plate solve succeeds.")
    print()

    image_height   = prompt_int(  "Image height           (pixels)",
                                  min_val=1)
    image_width    = prompt_int(  "Image width            (pixels)",
                                  min_val=1)
    horizon_pixel  = prompt_float("Horizon Y-pixel        (0=top of image)",
                                  min_val=0.0, max_val=float(image_height))
    fov_vertical   = prompt_float("Vertical FOV degrees   (from lens spec or astrometry.net)",
                                  min_val=0.001, max_val=180.0)
    print()
    print("  Rotation angle: degrees image +Y is rotated from North toward East.")
    print("  0 = north up. Auto-solve will override this if it runs.")
    rotation_angle = prompt_float("Field rotation angle   (degrees)",
                                  default=0.0, min_val=-360.0, max_val=360.0)

    dpp = fov_vertical / image_height
    return {
        "image_height":      image_height,
        "image_width":       image_width,
        "fov_vertical":      fov_vertical,
        "horizon_pixel":     horizon_pixel,
        "rotation_angle":    rotation_angle,
        "degrees_per_pixel": dpp,
    }


# ══════════════════════════════════════════════════════════════════════════════
#  STAR COLLECTION (manual)
# ══════════════════════════════════════════════════════════════════════════════

def _altitude_from_pixel(pixel_y, image_params):
    return (image_params["horizon_pixel"] - pixel_y) * image_params["degrees_per_pixel"]


def collect_stars_manual(image_params):
    print()
    print("  Enter matched stars from your plate solver.")
    print("  RA: decimal degrees 0-360.  Dec: decimal degrees -90 to +90.")
    print("  Pixel Y: Y-coordinate of the star (0=top). Min 2 stars; 4+ recommended.")
    print()

    n_stars = prompt_int("Number of stars", min_val=1, max_val=100)
    stars   = []

    for i in range(1, n_stars + 1):
        print(f"\n  -- Star {i} --")
        ra_deg  = prompt_float("  RA   (decimal degrees, 0-360)",    min_val=0.0,   max_val=360.0)
        dec_deg = prompt_float("  Dec  (decimal degrees, -90 to 90)",min_val=-90.0, max_val=90.0)
        pixel_y = prompt_float("  Y-pixel  (0=top)",
                               min_val=0.0, max_val=float(image_params["image_height"]))

        alt_apparent  = _altitude_from_pixel(pixel_y, image_params)
        refraction    = _refraction_correction(alt_apparent)
        alt_corrected = alt_apparent - refraction

        print(f"    Apparent alt : {alt_apparent:+.4f} deg")
        print(f"    Refraction   : -{refraction*60:.2f} arcmin")
        print(f"    Corrected alt: {alt_corrected:+.4f} deg")

        if alt_apparent < 0:
            print("    WARNING: below horizon -- check Y pixel.")
        elif alt_apparent < 5:
            print("    WARNING: below 5 deg -- refraction highly uncertain.")

        stars.append({
            "ra_deg":        ra_deg,
            "dec_deg":       dec_deg,
            "pixel_y":       pixel_y,
            "alt_apparent":  alt_apparent,
            "alt_corrected": alt_corrected,
            "source":        "manual",
        })

    return stars


# ══════════════════════════════════════════════════════════════════════════════
#  ASTRONOMICAL CALCULATIONS
# ══════════════════════════════════════════════════════════════════════════════

def julian_date(dt):
    """Meeus Chapter 7."""
    y = dt.year
    m = dt.month
    d = dt.day + (dt.hour + dt.minute / 60.0 + dt.second / 3600.0) / 24.0
    if m <= 2:
        y -= 1
        m += 12
    A = int(y / 100)
    B = 2 - A + int(A / 4)
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5


def gmst_degrees(jd):
    """Greenwich Mean Sidereal Time in degrees. Meeus Chapter 12."""
    T      = (jd - 2451545.0) / 36525.0
    gmst_s = (
        67310.54841
        + (876600.0 * 3600.0 + 8640184.812866) * T
        + 0.093104 * T ** 2
        - 6.2e-6   * T ** 3
    )
    return (gmst_s / 240.0) % 360.0


def _refraction_correction(alt_deg):
    """Bennett (1982) refraction in degrees. Subtract from apparent altitude."""
    if alt_deg < -1.0:
        return 0.0
    arcmin = 1.02 / math.tan(math.radians(alt_deg + 10.3 / (alt_deg + 5.11)))
    return max(0.0, arcmin / 60.0)


def predicted_altitude(lat_deg, lon_deg, ra_deg, dec_deg, gst_deg):
    """
    Geometric altitude via navigational triangle.
    sin(alt) = sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(LHA)
    LHA = GMST + observer_east_longitude - star_RA
    """
    lha = math.radians((gst_deg + lon_deg - ra_deg) % 360.0)
    lat = math.radians(lat_deg)
    dec = math.radians(dec_deg)
    sin_alt = (
        math.sin(lat) * math.sin(dec) +
        math.cos(lat) * math.cos(dec) * math.cos(lha)
    )
    return math.degrees(math.asin(max(-1.0, min(1.0, sin_alt))))


# ══════════════════════════════════════════════════════════════════════════════
#  GEOLOCATION SOLVER
# ══════════════════════════════════════════════════════════════════════════════

def _residuals(params, stars, gst_deg):
    lat, lon = params
    return sum(
        (predicted_altitude(lat, lon, s["ra_deg"], s["dec_deg"], gst_deg)
         - s["alt_corrected"]) ** 2
        for s in stars
    )


def solve_multistar(stars, gst_deg):
    """
    Multi-star least-squares via scipy L-BFGS-B with ten starting points.
    Returns (lat_deg, lon_deg, rms_deg).
    """
    starts = [
        ( 0.0,    0.0), ( 45.0,  -90.0), (-45.0,   90.0), ( 30.0, -100.0),
        ( 51.5,  -0.1), (-33.9,  151.2), ( 55.8,   37.6), ( 35.7,  139.7),
        ( 19.4,  -99.1), (-23.5, -46.6),
    ]
    best_result, best_fun = None, float("inf")
    for lat0, lon0 in starts:
        result = minimize(
            _residuals, x0=[lat0, lon0], args=(stars, gst_deg),
            method="L-BFGS-B",
            bounds=[(-90.0, 90.0), (-180.0, 180.0)],
            options={"ftol": 1e-14, "gtol": 1e-12, "maxiter": 20000}
        )
        if result.fun < best_fun:
            best_fun, best_result = result.fun, result
    lat, lon = best_result.x
    return lat, lon, math.sqrt(best_fun / len(stars))


def solve_single_star_fallback(star, gst_deg):
    """
    Single-star transit approximation.
    Returns (lat_south, lat_north, lon_approx).
    """
    dec = star["dec_deg"]
    alt = star["alt_corrected"]
    gha = (gst_deg - star["ra_deg"]) % 360.0
    lon = gha if gha <= 180.0 else gha - 360.0
    lat_s = max(-90.0, min(90.0, dec + (90.0 - alt)))
    lat_n = max(-90.0, min(90.0, dec - (90.0 - alt)))
    return lat_s, lat_n, lon


# ══════════════════════════════════════════════════════════════════════════════
#  OUTPUT
# ══════════════════════════════════════════════════════════════════════════════

def _fmt(deg, pos, neg):
    return f"{abs(deg):.5f} {pos if deg >= 0 else neg}  ({deg:+.5f})"


def print_result_multistar(lat, lon, rms, stars, gst_deg, method):
    section("Result -- Multi-Star Least-Squares Solution")
    print()
    print(f"  Plate solve      : {method}")
    print(f"  Latitude         : {_fmt(lat, 'N', 'S')}")
    print(f"  Longitude        : {_fmt(lon, 'E', 'W')}")
    print(f"  RMS residual     : {rms:.4f} deg  ({rms*60:.2f} arcmin)")
    print()
    header = "  #    Src       RA      Dec       Meas       Pred      Err     Err'"
    print(header)
    print("  " + "-" * 70)
    for i, s in enumerate(stars, 1):
        pred = predicted_altitude(lat, lon, s["ra_deg"], s["dec_deg"], gst_deg)
        err  = pred - s["alt_corrected"]
        src  = s.get("source", "")[:7]
        print(f"  {i:<4} {src:<8} {s['ra_deg']:>8.3f} {s['dec_deg']:>8.3f} "
              f"{s['alt_corrected']:>10.4f} {pred:>10.4f} {err:>8.4f} {err*60:>7.2f}")
    print()
    print(f"  Google Maps : https://www.google.com/maps?q={lat:.5f},{lon:.5f}")
    print(f"  Decimal     : {lat:.5f}, {lon:.5f}")
    print()

    if rms > 2.0:
        print("  HIGH RESIDUAL WARNING -- RMS > 2 degrees.")
        print("  Common causes: wrong horizon pixel, wrong observation time,")
        print("  plate solve on wrong field, stars below 10 deg.")
    elif rms > 0.5:
        print("  Moderate residual. Review stars with the largest individual errors.")
    else:
        print("  Residual is low. Solution is likely reliable.")


def print_result_single(lat_s, lat_n, lon, method):
    section("Result -- Single-Star Approximation  [LOW ACCURACY]")
    print()
    print(f"  Solve method : {method}")
    print()
    print("  WARNING: meridian transit assumed. Rarely exactly true.")
    print("  Install scipy for the multi-star optimizer: pip install scipy")
    print()
    print(f"  Solution A (star to your south) : {_fmt(lat_s, 'N', 'S')}")
    print(f"  Solution B (star to your north) : {_fmt(lat_n, 'N', 'S')}")
    print(f"  Longitude  (transit assumed)    : {_fmt(lon,   'E', 'W')}")
    print()
    print(f"  Maps A : https://www.google.com/maps?q={lat_s:.5f},{lon:.5f}")
    print(f"  Maps B : https://www.google.com/maps?q={lat_n:.5f},{lon:.5f}")


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    banner()

    # Step 1 -- image file
    section("Step 1 -- Image File")
    print()
    print("  Path to the photograph.")
    print("  Used for MAC timestamps and upload to the plate solver.")
    print("  Press Enter to skip (all values entered manually).")
    print()
    raw_path   = input("    Image file path: ").strip().strip('"').strip("'")
    image_path = raw_path or None

    if image_path and not os.path.exists(image_path):
        print(f"\n    [WARNING] File not found: {image_path}")
        print("    Continuing without the file.")
        image_path = None

    # Step 2 -- observation time
    obs_dt = collect_observation_datetime(image_path)
    print(f"\n    Using: {obs_dt.strftime('%Y-%m-%d  %H:%M:%S  UTC')}")

    # Step 3 -- image geometry (horizon pixel always required)
    image_params = collect_image_parameters()

    # Step 4 -- plate solve pipeline
    if image_path:
        stars, method = run_plate_solve(image_path, image_params)
    else:
        print()
        print("  No image file provided -- skipping automated plate solve.")
        section("Step 4 -- Star Data (Manual)")
        stars  = collect_stars_manual(image_params)
        method = SolveMethod.MANUAL

    # Step 5 -- compute
    section("Computing")
    jd  = julian_date(obs_dt)
    gst = gmst_degrees(jd)
    n   = len(stars)

    print(f"\n    Julian Date  : {jd:.6f}")
    print(f"    GMST         : {gst:.4f} deg  ({gst / 15.0:.4f} h)")
    print(f"    Stars used   : {n}")
    print(f"    Solve method : {method}")

    if SCIPY_AVAILABLE and n >= 2:
        print("\n    Running least-squares optimizer...")
        lat, lon, rms = solve_multistar(stars, gst)
        print_result_multistar(lat, lon, rms, stars, gst, method)
    else:
        if not SCIPY_AVAILABLE:
            print("\n    scipy unavailable -- using single-star fallback.")
        else:
            print("\n    Only one star -- using single-star fallback.")
        lat_s, lat_n, lon = solve_single_star_fallback(stars[0], gst)
        print_result_single(lat_s, lat_n, lon, method)

    section("Done")
    print()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n  Interrupted.")
        sys.exit(0)
