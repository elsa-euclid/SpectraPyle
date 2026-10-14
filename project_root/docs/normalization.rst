Normalization Methods
=====================

SpectraPyle distinguishes between **scalar normalization** and optional
**continuum normalization**. The six ``spectra_normalization`` strategies below
apply one multiplicative scale factor to each spectrum (or, for the Francis
method, derive relative scale factors from the ensemble). They preserve the
continuum shape of each individual spectrum.

Continuum normalization is a separate, optional preprocessing step controlled by
``continuum_normalize``. It divides each individual spectrum by a smooth fitted
continuum before resampling and stacking. This removes broad continuum-shape
differences and therefore changes the population weighting represented by the
stack.

Comparison Table
----------------

.. list-table::
   :header-rows: 1
   :widths: 20, 25, 30, 25

   * - Config Value
     - Function
     - Typical Use Case
     - Conserves
   * - ``"no_normalization"``
     - identity (with conservation mode)
     - Pre-normalized spectra or template stacking
     - Flux or luminosity
   * - ``"median"`` (default)
     - :func:`~spectraPyle.spectrum.normalization.normSpecMed`
     - General-purpose stacking
     - Flux
   * - ``"interval"``
     - :func:`~spectraPyle.spectrum.normalization.normSpecInterv`
     - Normalizing over a rest-frame wavelength range
     - Flux
   * - ``"integral"``
     - :func:`~spectraPyle.spectrum.normalization.normSpecIntegrMean`
     - Integrating over entire spectrum or region
     - Flux (integrated)
   * - ``"custom"``
     - :func:`~spectraPyle.spectrum.normalization.normSpecCustom`
     - Using external catalog values (e.g., photometry)
     - Flux (catalog-driven)
   * - ``"template"``
     - :func:`~spectraPyle.spectrum.normalization.francis1991_normalize`
     - Multi-spectrum normalization via iterative reference scaling
     - Flux (relative)

.. note::

   The scalar normalization functions return normalized flux, normalized error,
   and the scalar normalization value. Continuum normalization instead returns
   the wavelength-dependent fitted continuum in addition to normalized flux and
   error. The legacy ``norm_factor`` diagnostic continues to describe only the
   scalar normalization step.

No Normalization
----------------

**Config value:** ``"no_normalization"``

**Use case:** When spectra are already normalized externally, or when physical information must be preserved (e.g., measuring emission line fluxes) or when applying Francis et al. (1991) method (spectra are normalized following a dedicated algorithm).

This is a pure identity operation. Each spectrum stays unchanged. The ``conservation`` field must specify whether the user wants to **conserve the flux** (e.g., simpler but optimal for conserving emission line ratios) or to **conserve luminosity** (e.g., to preserve the intrinsic luminosity of the spectral features, for instance to estimate intrinsic properties from emission lines, such as star formation rates from the Hα line)  

.. code-block:: yaml

   norm:
     spectra_normalization: "no_normalization"
     conservation: "flux"

Median Normalization
--------------------

**Config value:** ``"median"`` (default)

**Function:** :func:`~spectraPyle.spectrum.normalization.normSpecMed`

**Formula:**

.. math::

   f_{\mathrm{norm}} = \frac{f}{\mathrm{median}(f)}

Each spectrum is divided by the median of its flux array:

.. math::

   \mathrm{median} = \mathrm{nanmedian}(f_{\lambda})

**Use case:** General-purpose stacking of spectra with widely varying flux levels. Median normalization is robust to outliers and provides a representative normalization for most spectroscopic surveys.

**YAML Example:**

.. code-block:: yaml

   norm:
     spectra_normalization: "median"

Interval Normalization
----------------------

**Config value:** ``"interval"``

**Function:** :func:`spectraPyle.spectrum.normalization.normSpecInterv`

**Definition**

.. math::

   f_{\mathrm{norm}} = \frac{f}{S}

The normalization factor :math:`S` is computed from the flux within a rest-frame wavelength interval
:math:`[\lambda_{\mathrm{min}}, \lambda_{\mathrm{max}}]`, using one of the following statistics:

.. math::

   S \in \{\mathrm{median}, \mathrm{mean}, \mathrm{maximum}, \mathrm{minimum}\}

**Description**

This function normalizes a spectrum over a user-defined rest-frame wavelength range by dividing the flux
by a representative statistic of that interval.

**Guidelines for choosing** :math:`S`

- :math:`\mathrm{median}` or :math:`\mathrm{mean}`  
  Recommended when the continuum is well-defined and free of strong emission or absorption features.

- :math:`\mathrm{maximum}`  
  Useful for normalizing to the peak of an emission line. When stacking spectra, this approach preserves
  the average full width at half maximum (FWHM) of the emission line.

- :math:`\mathrm{minimum}`  
  Analogous to ``maximum``, but suited for absorption features.


**YAML Example:**

.. code-block:: yaml

   norm:
     spectra_normalization: "interval"
     lambda_norm_min: 1400.0
     lambda_norm_max: 1600.0
     interval_norm_statistics: "median"

The allowed values for ``interval_norm_statistics`` are ``"median"``, ``"mean"``, ``"maximum"``, and ``"minimum"``.

Integral (Integrated Mean) Normalization
-----------------------------------------

**Config value:** ``"integral"``

**Function:** :func:`~spectraPyle.spectrum.normalization.normSpecIntegrMean`

**Formula:**

.. math::

   f_{\mathrm{norm}} = \frac{f}{\langle f \rangle_{\lambda}}

where the integrated mean flux is:

.. math::

   \langle f \rangle_{\lambda} = \frac{\sum_{i} f_i \Delta\lambda_i}{\lambda_{\mathrm{max}} - \lambda_{\mathrm{min}}}

This is the total integrated flux divided by the wavelength span, equivalent to averaging all flux values.

**Use case:** When you want a normalization that accounts for the entire spectrum's flux content rather than a single point or interval. Useful for broadband stacking where the average flux level is more meaningful than the median.

**YAML Example:**

.. code-block:: yaml

   norm:
     spectra_normalization: "integral"

Custom Normalization
---------------------

**Config value:** ``"custom"``

**Function:** :func:`~spectraPyle.spectrum.normalization.normSpecCustom`

**Formula:**

.. math::

   f_{\mathrm{norm}} = \frac{f}{N_{\mathrm{custom}}}

where :math:`N_{\mathrm{custom}}` is a per-spectrum scalar read from a catalog column.

**Use case:** When you have external reference fluxes (e.g., photometric measurements in a fixed band) that you want to use as the normalization reference. This is particularly useful for samples with heterogeneous data or when normalizing to a common photometric standard.

**YAML Example:**

.. code-block:: yaml

   norm:
     spectra_normalization: "custom"
   catalog_columns:
     custom_normalization:
       custom_column_name: "ref_flux_band_mag"

The catalog must include the column specified by ``custom_column_name``. Each value should be finite and positive.

Continuum Normalization (optional)
----------------------------------

**Config flag:** ``continuum_normalize: true``

**Function:** :func:`~spectraPyle.spectrum.normalization.normSpecContinuum`

Continuum normalization is independent of the scalar ``spectra_normalization``
choice. For every individual spectrum and every selected grism/slit, SpectraPyle
fits a smooth polynomial continuum :math:`C(\lambda)` and applies

.. math::

   f_{\mathrm{cont}}(\lambda) = \frac{f(\lambda)}{C(\lambda)}

and

.. math::

   \sigma_{\mathrm{cont}}(\lambda) = \frac{\sigma(\lambda)}{C(\lambda)}.

The operation is performed **before resampling**. It is implemented in the
generic spectrum-processing layer and therefore does not depend on Euclid, DESI,
or any other instrument module.

Why use it?
~~~~~~~~~~~

A scalar normalization, such as division by the median, does not change the
intrinsic equivalent width of an individual spectrum. However, spectra with
different continuum slopes can still contribute with different effective
weights at different wavelengths when they are stacked. Continuum normalization
removes most of this smooth shape dependence and produces a stack that is closer
to an average of the individual *fractional* spectral features.

This can be useful for absorption-line studies, template construction, or other
applications where continuum shape is not the quantity of interest. It should
not be enabled when the broad continuum or relative spectrophotometric shape is
part of the desired signal.

.. important::

   Continuum-normalized stacking does **not** guarantee that an equivalent width
   measured on the final stack is exactly equal to the arithmetic mean of the
   individual equivalent widths. Spectral coverage, masks, inverse-variance
   weighting, resolution, velocity broadening, and the continuum fit itself can
   all affect the statistic represented by the stack.

Fitting and feature protection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The continuum is fitted with an error-weighted polynomial. Iterative clipping is
performed on residuals normalized by the input 1-sigma errors, making the
clipping consistent with heteroscedastic spectra.

When ``continuum_mask_features: true`` (default), SpectraPyle also excludes:

- all absorption-feature bandpasses in the bundled absorption table; and
- the curated primary emission lines in the bundled vacuum line table.

The absorption pseudo-continuum sidebands are intentionally left available to
constrain the smooth continuum.

For stacks shifted to a common redshift, the feature mask is moved to the actual
stacking frame. This includes rest-frame, minimum-z, maximum-z, median-z, and
other non-zero common-frame stacks. In ``observed_frame`` mode SpectraPyle does
not require source redshifts, so it cannot place rest-frame features reliably;
automatic feature masking is therefore skipped and robust clipping is used
instead.

Polynomial order
~~~~~~~~~~~~~~~~

``continuum_norm_order`` can be either a single integer or a mapping keyed by
selected grism/slit name. The default is 2. Orders 1--2 are recommended as the
first choices because high-order polynomials can begin to absorb broad
astrophysical features. Higher orders should be validated for the specific
wavelength range and science case.

If the requested polynomial is insufficiently constrained, cannot be fitted, or
becomes non-positive over valid wavelengths, SpectraPyle automatically retries
with progressively lower orders, down to a positive constant median continuum.

**YAML examples:**

.. code-block:: yaml

   norm:
     spectra_normalization: median
     continuum_normalize: true
     continuum_norm_order: 2
     continuum_mask_features: true
     continuum_emission_mask_half_width: 12.0
     continuum_sigma_clip: 3.0

For an instrument with multiple spectral segments, the order may be set per
segment using the exact grism/slit names defined by that instrument module:

.. code-block:: yaml

   norm:
     spectra_normalization: median
     continuum_normalize: true
     continuum_norm_order:
       blue: 2
       red: 2
       default: 2

The mapping is instrument-agnostic: ``default``, ``merged``, ``arm1``, or any
other segment name can be used if that is how the corresponding instrument
module identifies its spectra.

Equivalent widths and Lick-style indices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A constant scalar normalization leaves an individual equivalent width exactly
unchanged. A smooth continuum division will usually preserve narrow local
line-to-continuum contrasts approximately, but it is not mathematically neutral
for every index definition. Lick-style indices use local pseudo-continua and are
often relatively insensitive to a low-order multiplicative correction; broad
molecular bands and wide indices are more sensitive and should be explicitly
validated.

Uncertainty caveat
~~~~~~~~~~~~~~~~~~

The returned error array is divided by the fitted continuum, but the uncertainty
of the continuum model itself is **not** propagated. This approximation is most
appropriate when the continuum is reasonably constrained. For very low
continuum-S/N spectra, compare results with and without continuum normalization
and inspect the stability of the fitted continuum.

Continuum normalization cannot currently be combined with
``spectra_normalization: template``. The Francis-style method derives
stack-level overlap scalings and represents a different normalization strategy.

Template Normalization (Francis et al. 1991)
---------------------------------------------

**Config value:** ``"template"``

**Function:** :func:`~spectraPyle.spectrum.normalization.francis1991_normalize`

**Method Overview:**

The Francis et al. (1991) method normalizes a set of spectra iteratively against a consensus reference spectrum. This is fundamentally different from per-spectrum normalization: instead of dividing each spectrum by a statistic computed from its own flux, each spectrum is scaled to match a shared reference.

**Algorithm:**

1. **Select an anchor spectrum** — the first spectrum with sufficient valid pixels (default: ≥ 50 pixels) is designated as the reference.

2. **Initialize the reference** — the anchor is normalized by its own median, and becomes the initial reference flux.

3. **Iterative scaling** — for each subsequent spectrum:

   - Compute the overlap region with the current reference (pixels with finite flux in both).
   - Estimate a scaling factor :math:`\alpha` by computing the median (or mean) of the ratio :math:`f_{\mathrm{ref}} / f_i` in the overlap region, with optional sigma-clipping.
   - Scale the spectrum: :math:`f_{\mathrm{norm}} = \alpha \cdot f_i`.
   - Update the reference using a count-weighted average: :math:`f_{\mathrm{ref,new}} = \frac{f_{\mathrm{ref}} \cdot c + f_{\mathrm{norm}}}{c + 1}` (where :math:`c` is the overlap count).
   - Add any new pixels from the current spectrum to the reference.

4. **Output** — normalized flux and error for all spectra.

**Use case:** Producing spectral templates of homogeneous sources. The iterative approach ensures that all spectra are mutually consistent while preserving the relative flux levels within overlapping regions. Also used when you need relative flux conservation across a heterogeneous sample.

**YAML Example:**

.. code-block:: yaml

   norm:
     spectra_normalization: "template"

.. seealso::

   - :func:`~spectraPyle.spectrum.normalization.francis1991_normalize` — full function documentation with parameters for sigma-clipping and overlap control.
   - **Francis, P. J., Hewett, P. C., Foltz, C. B., Chaffee, F. H., Weymann, R. J., & Morris, S. L.** (1991), *ApJ*, **373**, 465. Composite quasar spectra from the Hubble Space Telescope Faint Object Spectrograph.
