# stellar_geolocate.py

**Geolocate a photographer from a night-sky image using stellar astrometry.**

The tool determines where on Earth a photograph was taken by analyzing the visible star field. It plate-solves the image to extract celestial coordinates, then uses spherical astronomy and a multi-star least-squares optimizer to compute the observer's latitude and longitude. The result is a GPS coordinate and a direct Google Maps link.

---

## How It Works

Every night-sky photograph encodes the observer's position. The pattern of stars visible — which ones appear, at what altitude above the horizon, and in what orientation — is uniquely determined by four variables: latitude, longitude, date, and time. If the date and time are known (from file timestamps or manual entry), and the celestial coordinates of stars in the image can be identified (via plate solving), the observer's position can be solved as the latitude/longitude pair that minimizes the difference between predicted and measured stellar altitudes across all identified stars.

The core relationship is the navigational triangle altitude formula:

```
sin(alt) = sin(lat)·sin(dec) + cos(lat)·cos(dec)·cos(LHA)
```

where `LHA = GMST + observer_east_longitude − star_RA`. The tool solves for `lat` and `lon` simultaneously using a least-squares optimizer seeded from ten geographically distributed starting points to avoid local minima.

---

## Features

- Three-tier plate-solve pipeline with automatic fallback
- MAC timestamp extraction (Modified, Accessed, Created) with platform-aware behavior
- Bennett atmospheric refraction correction applied per star
- Multi-star least-squares geolocation via scipy L-BFGS-B
- Full gnomonic (TAN) WCS projection in pure Python — no FITS library required
- Optional astropy integration for full FITS WCS accuracy
- Per-star residual table with RMS quality assessment
- Single-star meridian-transit fallback when scipy is unavailable
- Google Maps link in output

---

## Requirements

**Python:** 3.8 or later

**Standard library only** is sufficient for manual-entry mode. Additional packages unlock automated pipeline stages:

| Package   | Required for                              | Install                  |
|-----------|-------------------------------------------|--------------------------|
| `requests`| nova API upload, WCS file download        | `pip install requests`   |
| `scipy`   | Multi-star least-squares optimization     | `pip install scipy`      |
| `astropy` | Full FITS WCS accuracy (optional upgrade) | `pip install astropy`    |

Install all at once:
```
pip install requests scipy astropy
```

The script detects which packages are available at startup and adjusts its behavior accordingly. It will not crash if optional packages are missing; it degrades gracefully and tells you exactly what is unavailable and why.

---

## Plate-Solve Pipeline

The tool attempts plate solving in the following order. Each tier is tried before falling back to the next.

### Tier 1 — nova.astrometry.net API

Requires `requests`. The script pings the API before asking anything; if the server is unreachable, this tier is silently skipped.

If reachable, the user is prompted for a free API key. The image is uploaded, and the script polls for job completion (up to 6 minutes, polling every 8 seconds). On success, calibration data and an optional WCS FITS file are downloaded. The FITS file is parsed with astropy if available; otherwise a gnomonic approximation is used. A 5×4 grid of sky positions is sampled from the WCS solution to build the star list automatically.

**Getting an API key:** Register a free account at [nova.astrometry.net](https://nova.astrometry.net), then find your key under *My Profile*.

### Tier 2 — Local solve-field

Requires the `solve-field` binary from the astrometry.net software package. The script checks PATH and several common installation directories. If found, the image is solved locally in a temporary directory. The resulting `.wcs` FITS file is parsed with astropy. This tier works entirely offline once the index files are installed.

**Installation:**
```bash
# Debian / Ubuntu
sudo apt install astrometry.net

# macOS (Homebrew)
brew install astrometry.net
```

Note: the astrometry.net package requires index files for blind solving. These are large downloads (several GB for full sky coverage). See [astrometry.net documentation](http://astrometry.net/doc/readme.html) for index file installation instructions.

### Tier 3 — Manual Entry

The fallback when both automated tiers fail or are declined. The user enters RA and Dec (decimal degrees) and pixel Y-coordinate for each star, sourced from any external plate solver such as [nova.astrometry.net](https://nova.astrometry.net) used manually, [KStars/EKOS](https://kstars.kde.org), [Siril](https://siril.org), or [AstroImageJ](https://www.astro.louisville.edu/software/astroimagej/).

A minimum of two stars is required. Four or more is strongly recommended; more stars produce lower residuals and more defensible results.

---

## Usage

```
python stellar_geolocate.py
```

The script is fully interactive. It walks through the following steps:

### Step 1 — Image File

Provide the path to the photograph. The file is used for two purposes: MAC timestamp extraction and upload to the plate solver. Press Enter to skip (all values will be entered manually).

### Step 2 — Observation Date and Time

If a file was provided, MAC timestamps are extracted and displayed. The user chooses which timestamp represents the capture time, or enters the date manually. All times are UTC.

**Platform behavior for file timestamps:**

| Platform | "Created" timestamp        | Notes                                  |
|----------|----------------------------|----------------------------------------|
| Windows  | `st_ctime` = creation time | Reliable creation time                 |
| macOS    | `st_birthtime`             | True birth time via HFS+/APFS metadata |
| Linux    | `st_ctime`                 | Inode change time, NOT creation time   |

On Linux, the Modified timestamp is generally the most useful for camera-transferred images. The script labels each timestamp clearly so there is no ambiguity.

**Accuracy note:** Time precision matters significantly. One minute of error in the observation time translates to approximately 0.25° of longitude error (~28 km at the equator). Use the most accurate timestamp available, or extract the original capture time from EXIF data before running the tool.

### Step 3 — Image Geometry

The following values are always required, regardless of solve method:

| Parameter          | Description                                                   |
|--------------------|---------------------------------------------------------------|
| Image height       | Total pixel height of the image                               |
| Image width        | Total pixel width of the image                                |
| Horizon Y-pixel    | Y-coordinate of the visible horizon line (0 = top of image)  |
| Vertical FOV       | Vertical field of view in degrees (from lens spec or solver)  |
| Field rotation     | Degrees image +Y is rotated from North toward East (default 0)|

The horizon pixel is the only value that cannot be determined by a plate solver — it is a scene geometry measurement that you must read by inspection from the image. In an image editor, place your cursor at the horizon line and note the Y-coordinate. If the plate solve succeeds, the FOV and rotation values are updated automatically from the WCS solution.

### Step 4 — Plate Solve

The three-tier pipeline runs automatically. If automated solving succeeds, a 5×4 grid of sky positions is sampled from the WCS solution; no star-by-star entry is needed. If falling back to manual entry, the user enters RA, Dec, and pixel Y for each star.

### Step 5 — Compute

Julian Date and Greenwich Mean Sidereal Time are computed from the observation datetime (Meeus, *Astronomical Algorithms*, Chapters 7 and 12). The multi-star optimizer runs and outputs a result.

---

## Output

```
── Result -- Multi-Star Least-Squares Solution ────────────────────────────────

  Plate solve      : nova.astrometry.net API
  Latitude         : 44.97832 N  (+44.97832)
  Longitude        : 93.26502 W  (-93.26502)
  RMS residual     : 0.0312 deg  (1.87 arcmin)

  #    Src       RA      Dec       Meas       Pred      Err     Err'
  ----------------------------------------------------------------------
  1    wcs_grid  83.822  -5.391    31.4821    31.4634  -0.0187   -1.12
  2    wcs_grid  95.988   9.221    38.2903    38.2711  -0.0192   -1.15
  ...

  Google Maps : https://www.google.com/maps?q=44.97832,-93.26502
  Decimal     : 44.97832, -93.26502
```

**RMS residual** is the root-mean-square difference between measured and predicted stellar altitudes across all stars. It is the primary quality indicator:

| RMS             | Interpretation                                              |
|-----------------|-------------------------------------------------------------|
| < 0.5°          | Low residual. Result is likely reliable.                    |
| 0.5° – 2.0°     | Moderate. Review stars with the largest individual errors.  |
| > 2.0°          | High. The result should not be trusted without investigation.|

Common causes of high residuals: incorrect horizon pixel, wrong observation time, plate solve converged on the wrong field, or stars below 10° altitude where refraction modeling degrades.

---

## Astronomical Methods

### Altitude Measurement

A star's altitude above the horizon is derived from its Y-pixel position in the image:

```
altitude = (horizon_pixel − star_pixel_y) × (FOV / image_height)
```

Positive altitude indicates the star is above the horizon. Pixel Y increases downward in image coordinates, so stars higher in the frame (smaller Y) correctly produce larger positive altitudes.

### Atmospheric Refraction

Bennett's (1982) refraction formula is applied to each star's apparent altitude before solving:

```
refraction (arcmin) = 1.02 / tan(alt + 10.3 / (alt + 5.11))
```

The corrected (geometric) altitude is then used in the optimizer. Refraction is significant below 15° altitude and approaches 35 arcminutes at the horizon. Stars below 5° are flagged as unreliable; stars below −1° receive no correction.

### WCS Projection

When a plate solve succeeds, pixel coordinates are converted to RA/Dec using the full gnomonic (TAN) projection. If astropy is available, the FITS WCS is parsed directly, capturing all distortion terms encoded by the solver. If not, the gnomonic transformation is computed from calibration parameters (plate scale, orientation, parity, field center) in pure Python. The two methods agree closely at field center and diverge slightly toward the image edges.

### Geolocation Optimizer

The least-squares objective function minimizes the sum of squared altitude residuals:

```
Σ ( predicted_altitude(lat, lon, RA_i, Dec_i, GMST) − measured_altitude_i )²
```

The optimizer is seeded from ten geographically distributed starting points and the best result (lowest residual) is returned. Bounds are enforced: latitude −90° to +90°, longitude −180° to +180°.

---

## Limitations and Caveats

**Time accuracy is critical.** One minute of error in UTC observation time produces approximately 0.25° (~28 km) of longitude error. Camera clocks drift; if the image was taken with a phone, the timestamp is likely accurate to within a second. Standalone cameras may be off by minutes or more if the clock was never set correctly.

**The horizon pixel must be measured manually.** No plate solver can determine where the horizon falls in the image. Errors in the horizon pixel propagate directly into altitude measurements for all stars and inflate the RMS residual uniformly.

**Single-star mode is inherently ambiguous.** Without scipy, the tool falls back to assuming the star is at meridian transit (LHA = 0), which is rarely true. Two latitude solutions are returned because the spherical geometry is symmetric. Treat single-star output as a rough search area, not a location.

**Refraction below 10° is poorly modeled.** Bennett's formula becomes increasingly uncertain at low altitudes. Avoid relying on stars below 10° if possible; the tool warns when stars fall below 5°.

**The tool establishes possible location, not identity.** A geographic coordinate is not by itself attribution. Time zone analysis, terrain matching, and corroborating metadata should be applied in any investigative context.

---

## Dependencies Summary

| Package    | Version   | Purpose                              | Required?     |
|------------|-----------|--------------------------------------|---------------|
| Python     | ≥ 3.8     | Runtime                              | Yes           |
| `requests` | any       | HTTP for nova API and WCS download   | For Tier 1/2  |
| `scipy`    | any       | L-BFGS-B least-squares optimizer     | Strongly recommended |
| `astropy`  | any       | Full FITS WCS parsing                | Optional      |

---

## Version History

| Version | Changes |
|---------|---------|
| 1.0     | Initial script; hardcoded values, single-star only, incorrect spherical trig |
| 2.0     | Corrected navigational triangle formula; multi-star scipy optimizer; Bennett refraction; MAC timestamp extraction; all values prompted interactively |
| 3.0     | Three-tier plate-solve pipeline (nova API → local solve-field → manual); full gnomonic WCS in pure Python; astropy FITS WCS integration; WCS grid sampling; automatic FOV update from plate scale |

---

## References

- Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed. Willmann-Bell.
- Bennett, G.G. (1982). The calculation of astronomical refraction in marine navigation. *Journal of Navigation*, 35(2), 255–259.
- Lang, D. et al. (2010). Astrometry.net: Blind astrometric calibration of arbitrary astronomical images. *The Astronomical Journal*, 139(5), 1782–1800.
