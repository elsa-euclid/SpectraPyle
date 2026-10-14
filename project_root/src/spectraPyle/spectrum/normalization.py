"""
Flux normalization for stacking.

Provides multiple normalization strategies controlled by ``config['norm_type']``,
including median-flux normalization over a wavelength range, and the
:func:`francis1991_normalize` template-based method (Francis et al. 1991).
"""

import re
from functools import lru_cache
from importlib import resources

import numpy as np
import pandas as pd

from spectraPyle.utils.log import get_logger
from spectraPyle.utils.exceptions import NormalizationError

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Feature masks for continuum normalization
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def _load_continuum_feature_bands():
    """Load bundled rest-frame spectral features used by the continuum mask.

    The absorption table contributes only the feature bandpasses (``Line limits``),
    while the pseudo-continuum sidebands remain available to constrain the fit.
    Emission-line wavelengths are masked symmetrically around their rest-frame
    centres. The tables are package data and are therefore available for every
    supported instrument, including the Generic FITS reader.

    Returns
    -------
    absorption_bands : tuple of tuple
        Rest-frame ``(lower, upper)`` absorption-feature limits in Angstrom.
    emission_centres : tuple of float
        Rest-frame emission-line centres in Angstrom.
    """

    absorption_bands = []
    emission_centres = []

    try:
        table = resources.files("spectraPyle").joinpath(
            "tables", "absorptions_table.csv"
        )
        with resources.as_file(table) as path:
            absorption = pd.read_csv(path)
        for value in absorption["Line limits"].dropna():
            numbers = re.findall(
                r"\d+(?:\.\d+)?",
                str(value).replace("–", "-").replace("—", "-"),
            )
            if len(numbers) >= 2:
                lo, hi = float(numbers[0]), float(numbers[1])
                absorption_bands.append((min(lo, hi), max(lo, hi)))
    except Exception as exc:  # pragma: no cover - package-data failure
        logger.warning(
            "Could not load absorption-feature table for continuum masking: %s",
            exc,
        )

    try:
        table = resources.files("spectraPyle").joinpath(
            "tables", "emission_lines_vacuum_table.csv"
        )
        with resources.as_file(table) as path:
            emission = pd.read_csv(path)
        # The bundled table contains many very weak transitions. Mask the
        # curated primary set (show == 1) by default; robust clipping handles
        # additional isolated lines without discarding most continuum pixels.
        if "show" in emission.columns:
            emission = emission[emission["show"].fillna(0).astype(int) == 1]
        emission_centres = [
            float(value)
            for value in emission["wavelength_AA"].dropna().to_numpy()
        ]
    except Exception as exc:  # pragma: no cover - package-data failure
        logger.warning(
            "Could not load emission-line table for continuum masking: %s",
            exc,
        )

    return tuple(absorption_bands), tuple(emission_centres)


def build_feature_mask(
    lbd,
    frame_scale=1.0,
    emission_halfwidth=12.0,
):
    """Build a mask of known spectral features for continuum fitting.

    Parameters
    ----------
    lbd : array-like
        Wavelength array in the frame in which the continuum is fitted.
    frame_scale : float, optional
        Multiplicative factor that maps bundled rest-frame feature wavelengths
        into ``lbd``. For a stack evaluated at redshift ``z_stack``, use
        ``1 + z_stack``. The default is 1 (rest frame).
    emission_halfwidth : float, optional
        Rest-frame half-width, in Angstrom, masked around every bundled
        emission line. The width is multiplied by ``frame_scale``.

    Returns
    -------
    ndarray of bool
        ``True`` for pixels excluded from the continuum fit.

    Notes
    -----
    The mask is instrument-independent: it uses only wavelength and bundled
    astrophysical line lists. In observed-frame stacking, where SpectraPyle does
    not require source redshifts, the processing pipeline intentionally does not
    apply this automatic rest-frame feature mask.
    """

    lbd = np.asarray(lbd, dtype=float)
    if not np.isfinite(frame_scale) or frame_scale <= 0:
        raise ValueError("frame_scale must be finite and > 0")
    if not np.isfinite(emission_halfwidth) or emission_halfwidth < 0:
        raise ValueError("emission_halfwidth must be finite and >= 0")

    absorption_bands, emission_centres = _load_continuum_feature_bands()
    mask = np.zeros(lbd.shape, dtype=bool)

    for lo, hi in absorption_bands:
        mask |= (lbd >= lo * frame_scale) & (lbd <= hi * frame_scale)

    halfwidth = emission_halfwidth * frame_scale
    for centre in emission_centres:
        mask |= np.abs(lbd - centre * frame_scale) <= halfwidth

    return mask


def _robust_scatter(values):
    """Return a robust 1-sigma scatter estimate for finite values."""

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan

    centre = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - centre))
    scatter = 1.4826 * mad
    if not np.isfinite(scatter) or scatter == 0:
        scatter = np.nanstd(values)
    return scatter


def normSpecContinuum(
    specid,
    lbd,
    flux,
    error,
    order=2,
    feature_mask=None,
    sigma_clip=3.0,
    max_iter=4,
    min_pix_per_param=5,
):
    """Continuum-normalize one spectrum with a robust polynomial fit.

    A low-order polynomial is fitted to continuum pixels using inverse-error
    weights. Known spectral features can be excluded with ``feature_mask`` and
    additional outliers are iteratively rejected in *error-normalized* residual
    space. Flux and 1-sigma error are then divided by the fitted continuum.

    If the requested polynomial is not sufficiently constrained, fails to fit,
    or becomes non-positive over valid wavelengths, progressively lower orders
    are attempted down to a constant positive median continuum.

    Parameters
    ----------
    specid : str or int
        Spectrum identifier, used only in log/error messages.
    lbd, flux, error : array-like
        Wavelength, flux, and 1-sigma uncertainty arrays for one spectrum.
    order : int, optional
        Requested polynomial order. Typical values are 1--2; higher orders
        should be validated for the science case.
    feature_mask : array-like of bool, optional
        Pixels to exclude from the fit (``True`` = excluded). If omitted,
        :func:`build_feature_mask` is evaluated in the rest frame.
    sigma_clip : float, optional
        Clipping threshold applied to robustly standardized residuals.
    max_iter : int, optional
        Maximum number of fit/clipping iterations per polynomial order.
    min_pix_per_param : int, optional
        Minimum number of usable continuum pixels per polynomial coefficient.

    Returns
    -------
    fluxNorm : ndarray
        Continuum-normalized flux.
    errNorm : ndarray
        Input 1-sigma error divided by the fitted continuum. Uncertainty in the
        fitted continuum itself is not included.
    continuum : ndarray
        Fitted continuum evaluated on ``lbd``.

    Raises
    ------
    ValueError
        If inputs are inconsistent or no positive continuum model can be found.
    """

    lbd = np.asarray(lbd, dtype=float)
    flux = np.asarray(flux, dtype=float)
    error = np.asarray(error, dtype=float)

    if not (lbd.shape == flux.shape == error.shape):
        raise ValueError(
            f"Continuum normalization failure for {specid}: wavelength, flux, "
            "and error arrays must have identical shapes."
        )
    if lbd.ndim != 1:
        raise ValueError(
            f"Continuum normalization failure for {specid}: expected 1-D arrays."
        )
    if not isinstance(order, (int, np.integer)) or isinstance(order, (bool, np.bool_)):
        raise ValueError("Continuum polynomial order must be an integer")
    if order < 0:
        raise ValueError("Continuum polynomial order must be >= 0")
    if not np.isfinite(sigma_clip) or sigma_clip <= 0:
        raise ValueError("sigma_clip must be finite and > 0")
    if max_iter < 1:
        raise ValueError("max_iter must be >= 1")
    if min_pix_per_param < 1:
        raise ValueError("min_pix_per_param must be >= 1")

    finite = (
        np.isfinite(lbd)
        & np.isfinite(flux)
        & np.isfinite(error)
        & (error > 0)
    )

    if feature_mask is None:
        feature_mask = build_feature_mask(lbd)
    feature_mask = np.asarray(feature_mask, dtype=bool)
    if feature_mask.shape != lbd.shape:
        raise ValueError(
            f"Continuum normalization failure for {specid}: feature_mask must "
            "match the spectrum shape."
        )

    continuum_pixels = finite & (~feature_mask)
    n_continuum = int(np.sum(continuum_pixels))
    if n_continuum < 2:
        raise NormalizationError(
            f"Continuum normalization failure for {specid}: too few valid "
            f"continuum pixels ({n_continuum})."
        )

    requested_order = int(order)
    last_reason = "continuum fit failed"

    for effective_order in range(requested_order, -1, -1):
        min_required = max(
            effective_order + 1,
            min_pix_per_param * (effective_order + 1),
        )
        if n_continuum < min_required and effective_order > 0:
            last_reason = (
                f"only {n_continuum} continuum pixels for order {effective_order}"
            )
            continue

        if effective_order == 0:
            continuum_value = np.nanmedian(flux[continuum_pixels])
            if np.isfinite(continuum_value) and continuum_value > 0:
                continuum = np.full(lbd.shape, continuum_value, dtype=float)
                if effective_order != requested_order:
                    logger.warning(
                        "Continuum normalization for %s: reduced polynomial "
                        "order %s -> 0 (%s).",
                        specid,
                        requested_order,
                        last_reason,
                    )
                return flux / continuum, error / continuum, continuum
            last_reason = f"invalid median continuum ({continuum_value})"
            continue

        fit_domain = continuum_pixels.copy()
        x0 = np.nanmedian(lbd[continuum_pixels])
        half_range = 0.5 * (
            np.nanmax(lbd[continuum_pixels]) - np.nanmin(lbd[continuum_pixels])
        )
        if not np.isfinite(half_range) or half_range <= 0:
            last_reason = "degenerate wavelength range"
            continue
        x = (lbd - x0) / half_range

        coefficients = None
        fit_failed = False
        for _ in range(max_iter):
            if np.sum(fit_domain) < min_required:
                fit_failed = True
                last_reason = "too few pixels after continuum clipping"
                break

            try:
                weights = 1.0 / error[fit_domain]
                coefficients = np.polyfit(
                    x[fit_domain],
                    flux[fit_domain],
                    effective_order,
                    w=weights,
                )
            except (np.linalg.LinAlgError, ValueError, FloatingPointError) as exc:
                fit_failed = True
                last_reason = f"polynomial fit failed: {exc}"
                break

            model = np.polyval(coefficients, x)
            pulls = (flux - model) / error
            centre = np.nanmedian(pulls[fit_domain])
            scatter = _robust_scatter(pulls[fit_domain])

            if not np.isfinite(scatter) or scatter == 0:
                break

            new_fit_domain = (
                continuum_pixels
                & np.isfinite(pulls)
                & (np.abs(pulls - centre) <= sigma_clip * scatter)
            )

            if np.array_equal(new_fit_domain, fit_domain):
                break
            if np.sum(new_fit_domain) < min_required:
                last_reason = "too few pixels after continuum clipping"
                break
            fit_domain = new_fit_domain

        if fit_failed or coefficients is None:
            continue

        # Refit once on the final accepted pixel set so the returned model is
        # exactly consistent with the final clipping mask.
        try:
            weights = 1.0 / error[fit_domain]
            coefficients = np.polyfit(
                x[fit_domain],
                flux[fit_domain],
                effective_order,
                w=weights,
            )
            continuum = np.polyval(coefficients, x)
        except (np.linalg.LinAlgError, ValueError, FloatingPointError) as exc:
            last_reason = f"final polynomial fit failed: {exc}"
            continue

        valid_model = finite & np.isfinite(continuum) & (continuum > 0)
        if np.all(valid_model[finite]):
            if effective_order != requested_order:
                logger.warning(
                    "Continuum normalization for %s: reduced polynomial order "
                    "%s -> %s (%s).",
                    specid,
                    requested_order,
                    effective_order,
                    last_reason,
                )
            return flux / continuum, error / continuum, continuum

        last_reason = "fitted continuum is non-positive or non-finite"

    raise NormalizationError(
        f"Continuum normalization failure for {specid}: no valid positive "
        f"continuum model found ({last_reason})."
    )

def normSpecInterv(specid, lbd, flux, error, lambdamin_norm, lambdamax_norm, norm_stat):
    """Normalize spectrum over a wavelength interval.

    Parameters
    ----------
    specid : str or int
        Spectrum identifier (for error messages).
    lbd : ndarray
        Wavelength array.
    flux : ndarray
        Flux array.
    error : ndarray
        Flux error array.
    lambdamin_norm : float
        Minimum wavelength of normalization interval (Å).
    lambdamax_norm : float
        Maximum wavelength of normalization interval (Å).
    norm_stat : {'median', 'mean', 'maximum', 'minimum'}
        Statistic to compute over the interval.

    Returns
    -------
    tuple
        - fluxNorm : ndarray
            Normalized flux
        - errNorm : ndarray
            Normalized error
        - norm : float
            Normalization value used

    Raises
    ------
    ValueError
        If normalization interval is outside wavelength range or contains no valid data.
    """

    vecnorm = flux[(lbd >= lambdamin_norm) & (lbd <= lambdamax_norm)]

    if (lambdamin_norm < np.nanmin(lbd)) or (lambdamax_norm > np.nanmax(lbd)):
        logger.warning(
            f"Spectrum {specid}: normalization interval "
            f"[{lambdamin_norm}, {lambdamax_norm}] Å partially outside the "
            f"spectrum wavelength range [{np.round(lbd[0], 2)}, "
            f"{np.round(lbd[-1], 2)}] Å."
        )

    if len(vecnorm) == 0:
        raise NormalizationError(
            f"normalization interval [{lambdamin_norm}, {lambdamax_norm}] Å "
            f"outside the spectrum wavelength range "
            f"[{np.round(lbd[0], 2)}, {np.round(lbd[-1], 2)}] Å"
        )

    if np.all(np.isnan(vecnorm)):
        raise NormalizationError(
            f"no valid flux values in normalization interval "
            f"[{lambdamin_norm}, {lambdamax_norm}] Å"
        )
    if norm_stat == 'median':
        norm = np.nanmedian(vecnorm)
    elif norm_stat == 'mean':
        norm = np.nanmean(vecnorm)
    elif norm_stat == 'maximum':
        norm = np.nanmax(vecnorm)
    elif norm_stat == 'minimum':
        norm = np.nanmin(vecnorm)
    else:
        raise ValueError(
            f"Normalization statistic '{norm_stat}' is not supported."
        )

    if (not np.isfinite(norm)) or (norm <= 0):
        raise NormalizationError(
            f"invalid normalization value ({norm}) in interval "
            f"[{lambdamin_norm}, {lambdamax_norm}] Å"
        )

    fluxNorm = flux / norm
    errNorm = error / norm
    return fluxNorm, errNorm, norm


def normSpecMed(lbd, flux, error):
    """Normalize spectrum by median flux.

    Parameters
    ----------
    lbd : ndarray
        Wavelength array (unused, for API consistency).
    flux : ndarray
        Flux array.
    error : ndarray
        Flux error array.

    Returns
    -------
    tuple
        - fluxNorm : ndarray
            Normalized flux
        - errNorm : ndarray
            Normalized error
        - norm : float
            Median flux value used for normalization
    """
    norm = np.nanmedian(flux)

    if (not np.isfinite(norm)) or (norm <= 0):
        raise NormalizationError(
            f"invalid median normalization value ({norm})"
        )

    fluxNorm = flux / norm
    errorNorm = error / norm
    return fluxNorm, errorNorm, norm

def normSpecIntegrMean(lbd, flux, error):
    """Normalize spectrum by integrated mean flux.

    Parameters
    ----------
    lbd : ndarray
        Wavelength array.
    flux : ndarray
        Flux array.
    error : ndarray
        Flux error array.

    Returns
    -------
    tuple
        - fluxNorm : ndarray
            Normalized flux
        - errNorm : ndarray
            Normalized error
        - norm : float
            Integrated mean flux value
    """
    diff = np.full_like(lbd, 0, dtype=float)
    diff[:-1] = np.diff(lbd)
    diff[-1] = diff[-2]
    norm = np.nansum(flux*diff) / (lbd[-1]-lbd[0]) ## option n.1
    #norm = np.nansum(flux) / len(lbd) ## option n.2 equivalent to option n1

    if (not np.isfinite(norm)) or (norm <= 0):
        raise NormalizationError(
            f"invalid integral normalization value ({norm})"
        )

    fluxNorm = flux / norm
    errorNorm = error / norm
    return fluxNorm, errorNorm, norm

def normSpecCustom(lbd, flux, error, norm):
    """Normalize a spectrum by a caller-supplied scalar value.

    The normalization scalar ``norm`` is typically read from a catalog column
    (e.g. a photometric flux in a reference band). Use
    ``spectra_normalization: "custom"`` in the configuration and set
    ``catalog_columns.custom_normalization.custom_column_name`` to the column name.

    Parameters
    ----------
    lbd : array-like
        Wavelength array (Å). Accepted for API consistency but not used.
    flux : ndarray
        Flux array.
    error : ndarray
        Error (1-sigma) array.
    norm : float
        Normalization scalar. Must be finite and > 0.

    Returns
    -------
    fluxNorm : ndarray
        Normalized flux array.
    errorNorm : ndarray
        Normalized error array.
    norm : float
        The normalization value applied.

    Raises
    ------
    ValueError
        If ``norm`` is not finite or is ≤ 0.
    """
    if (not np.isfinite(norm)) or (norm <= 0):
        raise NormalizationError(
            f"invalid custom normalization value ({norm})"
        )

    fluxNorm = flux / norm
    errorNorm = error / norm
    return fluxNorm, errorNorm, norm

def francis1991_normalize(
    stackArr,
    stackArrErr,
    min_overlap=50,
    norm_stat="median",
    sigma_clip=3.0,
    max_iter_clip=3,
    feature_mask=None,
    eps=1e-10,
    spectrum_ids=None,
    eligible_mask=None,
    return_status=False,
):
    """
    Robust Francis-style normalization with sigma-clipped overlap scaling.

    Parameters
    ----------
    stackArr : ndarray (Npix, Nspec)
        Flux array
    stackArrErr : ndarray (Npix, Nspec)
        1-sigma errors
    min_overlap : int
        Minimum overlap for normalization
    norm_stat : {'median', 'mean'}
        Statistic used after clipping
    feature_mask : ndarray (Npix,) or None
        True = exclude from normalization
    sigma_clip : float
        Sigma clipping threshold (default: 3)
    max_iter_clip : int
        Max iterations for clipping
    eps : float
        Small value to avoid division issues

    Returns
    -------
    norm_flux : ndarray
    norm_err : ndarray
    alphas : ndarray (Nspec,) — per-spectrum normalization factors
    template_status : ndarray (Nspec,), optional
        Returned only when ``return_status=True``. Values identify spectra
        rejected by the template-normalization stage.
    """

    logger.info("Starting Francis 1991-like normalization")

    stat_func = np.nanmedian if norm_stat == "median" else np.nanmean

    Npix, Nspec = stackArr.shape

    norm_flux = np.full_like(stackArr, np.nan)
    norm_err = np.full_like(stackArrErr, np.nan)
    alphas = np.full(Nspec, np.nan)
    template_status = np.full(Nspec, "OK", dtype=object)

    if spectrum_ids is None:
        spectrum_ids = np.arange(Nspec)
    if eligible_mask is None:
        eligible_mask = np.ones(Nspec, dtype=bool)
    else:
        eligible_mask = np.asarray(eligible_mask, dtype=bool)
        if eligible_mask.shape != (Nspec,):
            raise ValueError("eligible_mask must have shape (Nspec,)")
    template_status[~eligible_mask] = "SKIP"

    def spectrum_label(i):
        value = spectrum_ids[i]
        if isinstance(value, bytes):
            value = value.decode("utf-8", errors="replace")
        return str(value)

    # --------------------------------------------------
    # Feature mask
    # --------------------------------------------------
    if feature_mask is None:
        feature_mask = np.zeros(Npix, dtype=bool)

    usable = (
        np.isfinite(stackArr)
        & np.isfinite(stackArrErr)
        & (~feature_mask[:, None])
    )

    counts = usable.sum(axis=0)

    # --------------------------------------------------
    # anchor = first (lowest-z) spectrum with enough pixels
    # --------------------------------------------------
    anchor = None
    alpha0 = np.nan
    for i in range(Nspec):
        if not eligible_mask[i]:
            continue
        if counts[i] < min_overlap:
            template_status[i] = "TEMPLATE_NO_OVERLAP"
            continue

        candidate = stat_func(stackArr[:, i][usable[:, i]])
        if not np.isfinite(candidate) or candidate <= 0:
            template_status[i] = "TEMPLATE_INVALID_NORMALIZATION"
            logger.warning(
                f"Spectrum {spectrum_label(i)} rejected during Francis 1991-like "
                f"normalization: invalid anchor normalization ({candidate})."
            )
            continue

        anchor = i
        alpha0 = candidate
        template_status[i] = "OK"
        break

    if anchor is None:
        raise ValueError(
            "Francis 1991-like normalization could not find a valid anchor spectrum"
        )

    logger.info(
        f"Francis 1991-like normalization: anchor spectrum "
        f"{spectrum_label(anchor)} (index {anchor})"
    )

    # --------------------------------------------------
    # Initialize reference
    # --------------------------------------------------
    f0 = stackArr[:, anchor]
    e0 = stackArrErr[:, anchor]

    valid0 = usable[:, anchor]

    norm_flux[:, anchor] = f0 / alpha0
    norm_err[:, anchor] = e0 / alpha0
    alphas[anchor] = alpha0

    ref_flux = norm_flux[:, anchor].copy()

    ref_count = np.zeros(Npix)
    ref_count[np.isfinite(ref_flux)] = 1

    # --------------------------------------------------
    # Robust alpha estimator
    # --------------------------------------------------

    def sig_clip_alpha(vals):
        for _ in range(max_iter_clip):
            mu = np.nanmedian(vals)
            std = np.nanstd(vals)

            if std == 0 or not np.isfinite(std):
                break

            good = np.abs(vals - mu) < sigma_clip * std

            if good.sum() == len(vals):
                break

            vals = vals[good]

            if vals.size < min_overlap:
                return np.nan
        return vals

    def robust_alpha(ref_vals, cur_vals):

        mask = (
            np.isfinite(ref_vals)
            & np.isfinite(cur_vals)
            & (np.abs(cur_vals) > eps)
        )

        if mask.sum() < min_overlap:
            return np.nan

        ## median of ratio
        #ratio = ref_vals[mask] / cur_vals[mask]
        #vals = ratio.copy()
        #vals = sig_clip_alpha(vals)
        #return stat_func(vals)

        ## ratio of medians
        ref_vals_sc = sig_clip_alpha(ref_vals[mask])
        cur_vals_sc = sig_clip_alpha(cur_vals[mask])

        ## with sigma clipping
        #return stat_func(ref_vals_sc) / stat_func(cur_vals_sc)

        ## no sigma clipping
        return stat_func(ref_vals[mask]) / stat_func(cur_vals[mask])

    # --------------------------------------------------
    # Incremental normalization
    # --------------------------------------------------
    for i in range(anchor + 1, Nspec):

        if not eligible_mask[i]:
            continue

        fi = stackArr[:, i]
        ei = stackArrErr[:, i]

        valid_i = np.isfinite(fi) & np.isfinite(ei)
        valid_r = np.isfinite(ref_flux)

        overlap = valid_i & valid_r & (~feature_mask)

        if overlap.sum() < min_overlap:
            template_status[i] = "TEMPLATE_NO_OVERLAP"
            logger.warning(
                f"Spectrum {spectrum_label(i)} rejected during Francis 1991-like "
                f"normalization: only {overlap.sum()} overlapping valid pixels "
                f"(minimum required: {min_overlap})."
            )
            continue

        ref_vals = ref_flux[overlap]
        cur_vals = fi[overlap]

        alpha = robust_alpha(ref_vals, cur_vals)

        if not np.isfinite(alpha) or (alpha <= 0) or (alpha > 1000):
            template_status[i] = "TEMPLATE_INVALID_NORMALIZATION"
            logger.warning(
                f"Spectrum {spectrum_label(i)} rejected during Francis 1991-like "
                f"normalization: invalid scale factor ({alpha})."
            )
            alphas[i] = np.nan
            norm_flux[:, i] = np.nan
            norm_err[:, i] = np.nan
            continue

        # normalize
        norm_flux[:, i] = alpha * fi
        norm_err[:, i] = alpha * ei
        alphas[i] = alpha

        # update reference (count-weighted)
        both = overlap

        ref_flux[both] = (
            ref_flux[both] * ref_count[both] + norm_flux[both, i]
        ) / (ref_count[both] + 1)

        ref_count[both] += 1

        # add new pixels
        new_only = valid_i & (~valid_r)

        ref_flux[new_only] = norm_flux[new_only, i]
        ref_count[new_only] = 1

    if return_status:
        return norm_flux, norm_err, alphas, np.asarray(template_status, dtype=str)
    return norm_flux, norm_err, alphas
