# SpectraPyle

**Stack thousands of astronomical spectra, uncover hidden signals.**

![Python](https://img.shields.io/badge/python-3.9%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![PyPI](https://img.shields.io/pypi/v/spectraPyle)

SpectraPyle is a flexible, scalable tool for stacking galaxy spectra. It shifts spectra
to a common rest frame, normalizes, resamples, sigma-clips, and combines them using median,
mean, weighted mean, or geometric mean statistics. Output is a FITS file ready for
scientific analysis.

**Supported instruments:** Euclid (NISP, grisms `red`/`blue`) · DESI (`merged`, back-compatibility only) · Generic (any standard FITS, grism `default`)

**Documentation:** Full documentation, tutorials, configuration reference, and API documentation are available on [Read the Docs](https://spectrapyle.readthedocs.io/).

---

## Quick Start

### Install

Install the latest stable release from PyPI:

```bash
pip install spectraPyle
```

This installs SpectraPyle together with the dependencies required by its interactive GUI.

If you also need a local JupyterLab environment and Voilà, install the notebook extras:

```bash
pip install "spectraPyle[notebook]"
```

For development or an editable installation from the source repository, see the [Development](#development) section below.

### Try SpectraPyle immediately

SpectraPyle includes a ready-to-use sample of 100 Euclid NISP Q1 spectra together with a pre-configured `default.gui`. This quick-start sample is intended for users who want to see the complete stacking workflow and its output immediately, before configuring their own data.

In a Jupyter notebook, launch the GUI with:

```python
from spectraPyle.gui import start

start()
```

When SpectraPyle is launched for the first time in a working directory, it creates:

```text
configs/GUI/default.gui
```

Load this configuration from the GUI and run SpectraPyle. The bundled catalogue and example spectra are located automatically: no input paths need to be configured.

The example performs a complete stack of the supplied Euclid spectra and generates the standard SpectraPyle outputs, including the interactive stacked-spectrum visualization.

The example data are installed as read-only package resources and are not copied into your project directory. Output files are instead written to the local `./outputs` directory.

Once you have verified the workflow, you can modify the GUI configuration for your own catalogue and spectra and save project-specific GUI, JSON, or YAML configurations under `./configs/`.

### Run via GUI

The GUI can be launched directly from any Jupyter environment:

```python
from spectraPyle.gui import start

start()
```

It is rendered inline in the notebook using `ipywidgets`; no local server is required.

On **ESA DataLabs**, SpectraPyle can therefore be used directly from a notebook. For Euclid internal data, see the [DataLabs documentation](https://spectrapyle.readthedocs.io/) for the dedicated workflow and configuration details.

For a standalone local Voilà interface, the source repository also provides:

```bash
python notebooks/run_gui.py
```

when SpectraPyle has been installed from source with the notebook extras.

### Run via CLI

After installation, use the command-line entry point:

```bash
spectrapyle --config path/to/config.yaml
```

Override individual config keys at runtime:

```bash
spectrapyle --config config.yaml --instrument.grisms '["red","blue"]'
```

### Minimal YAML config example

```yaml
instrument:
  instrument_name: euclid
  survey_name: wide
  grisms:
    - red
  data_release: Q1

io:
  input_dir: /path/to/catalog
  filename_in: my_catalog
  filename_in_extention: fits
  output_dir: /path/to/output
  grism_io:
    red:
      spectra_dir: /path/to/spectra

redshift:
  z_type: rest_frame

norm:
  norm_type: interval
  norm_range: [5100, 5500]

resampling:
  pixel_resampling_type: lambda
  pixel_size_type: instrumental
  nyquist_sampling: 5.0

sigmaclip:
  sigma_clipping_conditions: 3.0
```

---

## Configuration

Configs are YAML or JSON. All inputs pass through a strict validation pipeline:

```text
Widgets / JSON / YAML / CLI
    ↓
normalize_raw_config()        runtime/runtime_adapter.py  — empty strings→None, path coercion
    ↓
StackingConfig (Pydantic v2)  schema/schema.py             — typed field-level validation
    ↓
StackingConfigResolver        schema/schema.py             — cross-field rules, instrument defaults
    ↓
flatten_schema_model()        runtime/runtime_adapter.py    — back to flat dict (legacy shim)
    ↓
Stacking(flat_dict).run()     stacking/stacking.py
```

### Key config blocks

| Block | Purpose |
|---|---|
| `instrument` | `instrument_name` (`euclid`/`desi`/`generic`), `grisms` list |
| `io` / `grism_io` | Input catalog, output dir, per-grism spectra paths |
| `redshift` | `z_type` (`rest_frame`, `observed_frame`, `custom`, …), `z_value` |
| `norm` | Scalar normalization plus optional continuum normalization |
| `resampling` | Pixel grid type and size |
| `sigmaclip` | Sigma-clipping thresholds |
| `bootstrap` | Bootstrap uncertainty estimation |
| `parallel` | Multiprocessing CPU fraction |
| `plot` | Plotly output options |

### Optional continuum normalization

Continuum normalization is a separate, instrument-independent preprocessing
step. It fits each individual spectral segment before resampling and divides the
flux and error by the fitted smooth continuum.

```yaml
norm:
  spectra_normalization: median
  continuum_normalize: true
  continuum_norm_order: 2
  continuum_mask_features: true
  continuum_emission_mask_half_width: 12.0
  continuum_sigma_clip: 3.0
```

The scalar normalization does not change the intrinsic EW of an individual
spectrum. Continuum normalization instead removes smooth spectral-shape
differences before stacking; this can be useful for absorption-line work but
does not imply that an EW measured on the final stack is exactly the mean of
individual EWs. See the full documentation for details and caveats.

### Multi-grism

```yaml
instrument:
  grisms: ["red", "blue"]

io:
  grism_io:
    red:
      spectra_dir: /path/to/red_spectra
    blue:
      spectra_dir: /path/to/blue_spectra
```

Old single-grism configs (`grism_type: str`) are auto-migrated by `migrate_v1_to_v2()` in
`runtime/runtime_adapter.py`.

---

## Output Format

SpectraPyle writes a `.fits` file to `output_dir`.

### HDU[0] — metadata header

```python
from astropy.io import fits

with fits.open("result.fits") as hdu:
    redshift = hdu[0].header["REDSHIFT"]
```

### HDU[1] — stacked spectra table

```python
with fits.open("result.fits") as hdu:
    data = hdu[1].data
    wave = data.field("wavelength")
    spec = data.field("specMean")
```

| Column | Description |
|---|---|
| `wavelength` | Wavelength array (Å) |
| `specMean` / `specMeanDispersion` / `specMeanError` | Arithmetic mean stack + 1σ dispersion + 1σ uncertainty |
| `specMedian` / `specMedianDispersion` / `specMedianError` | Median stack (robust dispersion via MAD × 1.4826) |
| `specGeometricMean` / `specGeometricMeanDispersion` / `specGeometricMeanError` | Geometric mean (lenient) — non-positive flux values are excluded per pixel; result is NaN if all pixels have non-positive flux |
| `specStrictGeometricMean` / `specStrictGeometricMeanDispersion` / `specStrictGeometricMeanError` | Geometric mean (strict) — returns NaN for any pixel with ≥1 non-positive flux value; use when strict positivity is required |
| `specMode` / `specModeDispersion` / `specModeError` | Mode (Half-Sample Mode estimator, parameter-free); dispersion via MAD × 1.4826 |
| `specWeightedMean` / `specWeightedMeanDispersion` / `specWeightedMeanError` | Inverse-variance weighted mean |
| `spec16th`, `spec84th`, `spec98th`, `spec99th` | 16th, 84th, 97.73rd (≈2σ), and 99.73rd (≈3σ) percentiles |
| `initialPixelCount` | Total spectra submitted to the stack |
| `goodPixelCount` | Good flux values used per pixel |
| `badPixelCount` | Bad (masked) pixels per pixel |
| `sigmaClippedCount` | Values removed by sigma-clipping per pixel |
| `geomMeanPixelCount` | Spectra with positive flux that contributed to geometric mean per pixel |
| `templateNormMaskedCount` | Pixels masked during template normalization |

---

## Visualization

### Stacked spectrum plot

After a run, SpectraPyle automatically generates an interactive Plotly figure with two panels: pixel counts (top) and all six stacked estimators with error/dispersion bands, percentile envelopes, and astrophysical line markers (bottom).

```python
from spectraPyle.plot.plot import plotting

plotting("result_STACKING.fits", width=1200, height=700)
```

### 2D flux distribution — H5 array viewer

The intermediate HDF5 file (`*_array.h5`) stores the full matrix of individual resampled spectra before combination.

Use `plot_h5_heatmap()` to inspect the spectral distribution after a run:

```python
from spectraPyle.plot.plot import plot_h5_heatmap

fig = plot_h5_heatmap(
    h5_path="result_array.h5",
    fits_path="result_STACKING.fits",
    metric="specMedian",
    mode="heatmap",   # or "lines"
)

fig.show()
```

**Binning:** by default, one bin per original wavelength pixel with pixel-exact half-pixel edges. This preserves the physical spacing of linear and log-linear grids without any flux redistribution. Reducing `nbinsx` below `n_pixels` uses uniform bins and prints a warning.

---

## Supported Instruments

| Instrument | Grisms | Spectra modes |
|---|---|---|
| Euclid (NISP) | `red`, `blue` | `individual fits`, `combined fits`, `metadata path` |
| DESI *(back-compat)* | `merged` | `individual fits`, `combined fits`, `metadata path` |
| Generic | `default` | `individual fits`, `combined fits` |

> **DESI note:** The DESI driver is maintained for back-compatibility with Euclid–DESI crossmatched samples. It is not under active development. For new DESI-only workflows, consider using the Generic driver.

> **ESA DataLabs users:** If you are working on the DataLabs platform with Euclid internal data, use `spectra_mode: metadata path` and map the three SIR catalogue columns (`datalabs_path`, `file_name`, `hdu_index`) to SpectraPyle's metadata configuration. Access to internal data releases requires DataLabs authentication. See the full documentation on [Read the Docs](https://spectrapyle.readthedocs.io/) for the complete workflow.

---

## Architecture

### Config pipeline

```text
Widgets / JSON / YAML / CLI
    ↓
normalize_raw_config()        runtime/runtime_adapter.py
    ↓
StackingConfig (Pydantic v2)  schema/schema.py
    ↓
StackingConfigResolver        schema/schema.py
    ↓
flatten_schema_model()        runtime/runtime_adapter.py  ← legacy shim, to be deprecated
    ↓
Stacking(flat_dict).run()     stacking/stacking.py
```

### Key modules

| Module | Responsibility |
|---|---|
| `schema/schema.py` | Pydantic models: `StackingConfig`, `StackingConfigResolver`. Instrument rules from `instruments_rules.json`. |
| `runtime/runtime_adapter.py` | Input normalization, GUI→schema adapter, JSON/YAML loaders, version migration. |
| `stacking/stacking.py` | `Stacking` class: orchestrates the full run. Chunks spectra, writes HDF5, sigma-clips, computes statistics, saves FITS. |
| `process/processes.py` | `main_parallel()`: multiprocessing loop — read → shift → resample → normalize. |
| `instruments/euclid.py` | Euclid spectral reader and instrument-specific handling. |
| `instruments/desi.py` | DESI driver (back-compatibility with Euclid crossmatched samples only). |
| `instruments/generic.py` | Generic FITS driver: binary table + WCS header fallback, auto-detects wavelength range. |
| `instruments/instruments_rules.json` | Authoritative per-instrument constants, quality defaults, grisms, data releases. |
| `io/IO.py` | Catalog reader, wavelength grid builder, FITS output writer. |
| `io/filename_builder.py` | Auto-generates `filename_out` from config fields. |
| `spectrum/spectra.py` | Per-spectrum read, redshift shift, extinction correction. |
| `spectrum/resampling.py` | Wavelength grid resampling. |
| `spectrum/normalization.py` | Flux normalization; `francis1991_normalize` for template mode. |
| `statistic/statistics.py` | Stacking statistics and bootstrap uncertainties. |
| `plot/plot.py` | Plotly interactive visualization of the final stack and HDF5 spectral distribution. |
| `physics/extinction.py` | Galactic extinction correction (`dust_extinction`). |
| `utils/log.py` | Unified logging. |

### Adding a new instrument

1. Add an entry to `instruments/instruments_rules.json` with constants, quality defaults, surveys, grisms, and data releases.
2. Create `instruments/<name>.py` implementing the required spectral-reading and instrument-specific functions.
3. Register any instrument- or release-specific filename conventions.

---

## Development

Clone the repository and install SpectraPyle in editable mode with all optional dependencies:

```bash
git clone https://github.com/squai/SpectraPyle.git
cd SpectraPyle/project_root
pip install -e ".[all]"
```

### Lint and format

```bash
ruff check src/
black src/
```

### Sanity check

```bash
spectrapyle --config path/to/config.yaml
```

### Generate documentation locally

The latest documentation is built automatically on Read the Docs after updates to the repository. To build and preview it locally:

```bash
cd docs
make html
```

The generated documentation is available under `docs/_build/html/`.

---

## Citing SpectraPyle

If you use SpectraPyle in your research, please cite:

> Euclid Collaboration, Quai S., Pozzetti L., et al. 2026, *Astronomy & Astrophysics*, 707, A232.  
> DOI: [10.1051/0004-6361/202557329](https://doi.org/10.1051/0004-6361/202557329)

**BibTeX:**

```bibtex
@ARTICLE{Quai2026,
  author = {{Euclid Collaboration: Quai}, S. and {Pozzetti}, L. and {Talia}, M. and others},
  title = "{Euclid preparation. LXXXII. Predicting star-forming galaxy scaling relations with the spectral stacking code SpectraPyle}",
  journal = {\aap},
  year = 2026,
  month = mar,
  volume = {707},
  pages = {A232},
  doi = {10.1051/0004-6361/202557329},
  archivePrefix = {arXiv},
  eprint = {2509.16120},
  primaryClass = {astro-ph.GA}
}
```

Publications using SpectraPyle **must** also include the following acknowledgment:

> ELSA: Euclid Legacy Science Advanced analysis tools (Grant Agreement no. 101135203) is funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or Innovate UK. Neither the European Union nor the granting authority can be held responsible for them. UK participation is funded through the UK Horizon guarantee scheme under Innovate UK grant 10093177.

If you use SpectraPyle with **Euclid** data, you must also follow the [Euclid Consortium publication policy](https://www.euclid-ec.org/) and include the standard Euclid acknowledgment text.

For the full reference and detailed acknowledgment requirements, see the [SpectraPyle documentation](https://spectrapyle.readthedocs.io/).

---

## Contributors & Development Team

**Lead Developer & Maintainer:**
- **Salvatore Quai** (University of Bologna) — salvatore.quai@unibo.it

**Contributors:**
- Lucia Pozzetti
- Margherita Talia
- Zhiying Mao
- Xavier Lopez Lopez
- Marcel Lindmeir
- Elisabeta Lusso
- Sotiria Fotopoulou
- Michele Moresco

For inquiries or contributions, please contact Salvatore Quai.

---

## Contributing & Collaboration

We welcome discussions about new features, extensions, and scientific collaborations.

If you have ideas for:
- **New stacking statistics or methods** — we'd like to hear about them
- **Support for additional instruments or data formats** — let's discuss
- **Performance improvements or refactoring** — contributions are valued
- **Joint research projects using SpectraPyle** — we're open to collaboration

**Please open an issue or pull request** on the [GitHub repository](https://github.com/squai/SpectraPyle), or contact the development team:
- **salvatore.quai@unibo.it**, **salvatore.quai@gmail.com** (lead developer)

All contributions are reviewed and credited. We follow standard open-source practices for code review and testing.

---

## License

MIT — see `LICENSE`.

## Contact

Salvatore Quai — salvatore.quai@unibo.it, salvatore.quai@gmail.com
