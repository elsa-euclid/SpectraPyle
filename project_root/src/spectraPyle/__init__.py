"""
SpectraPyle — spectral stacking for Euclid and DESI.

Subpackages
-----------
stacking
    Full pipeline orchestrator (:class:`~spectraPyle.stacking.stacking.Stacking`).
schema
    Pydantic validation models.
runtime
    Config normalization and adapter layer.
instruments
    Per-instrument drivers (Euclid, DESI).
io
    Catalog I/O and FITS output.
spectrum
    Per-spectrum processing (shift, resample, normalize).
statistic
    Stacking estimators and bootstrap.
plot
    Plotly visualization.
physics
    Galactic extinction correction.
utils
    Logging utilities.
"""

__author__ = "Salvatore Quai"

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("spectraPyle")
except PackageNotFoundError:
    __version__ = "unknown"
