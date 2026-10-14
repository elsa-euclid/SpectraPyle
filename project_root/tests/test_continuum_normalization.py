import numpy as np
import pytest

from spectraPyle.spectrum.normalization import (
    _load_continuum_feature_bands,
    build_feature_mask,
    normSpecContinuum,
)


def test_bundled_feature_tables_are_available():
    absorption_bands, emission_centres = _load_continuum_feature_bands()
    assert len(absorption_bands) > 0
    assert len(emission_centres) > 0


def test_feature_mask_follows_stacking_frame():
    halpha = 6564.632
    z_stack = 0.5
    shifted = halpha * (1.0 + z_stack)
    wave = np.linspace(shifted - 100.0, shifted + 100.0, 401)

    mask = build_feature_mask(wave, frame_scale=1.0 + z_stack, emission_halfwidth=12.0)

    assert mask[np.argmin(np.abs(wave - shifted))]
    assert not mask[0]
    assert not mask[-1]


def test_continuum_normalization_recovers_quadratic_shape():
    wave = np.linspace(4000.0, 7000.0, 3001)
    x = (wave - 5500.0) / 1500.0
    true_continuum = 2.0 + 0.25 * x + 0.12 * x**2
    flux = true_continuum.copy()
    error = np.full_like(wave, 0.02)

    # Protect one absorption-like feature explicitly and let clipping reject
    # an isolated unmasked outlier elsewhere.
    feature_mask = (wave > 5150.0) & (wave < 5200.0)
    flux[feature_mask] *= 0.7
    flux[np.argmin(np.abs(wave - 6200.0))] += 2.0

    flux_norm, err_norm, fitted = normSpecContinuum(
        "synthetic",
        wave,
        flux,
        error,
        order=2,
        feature_mask=feature_mask,
    )

    continuum_pixels = ~feature_mask
    assert np.nanmedian(flux_norm[continuum_pixels]) == pytest.approx(1.0, abs=1e-6)
    assert np.nanmax(
        np.abs(fitted[continuum_pixels] / true_continuum[continuum_pixels] - 1.0)
    ) < 1e-6
    np.testing.assert_allclose(err_norm, error / fitted)


def test_requested_order_falls_back_when_underconstrained():
    wave = np.linspace(5000.0, 5014.0, 8)
    flux = np.linspace(1.0, 1.2, 8)
    error = np.full_like(wave, 0.05)

    _, _, fitted = normSpecContinuum(
        "few-pixels",
        wave,
        flux,
        error,
        order=2,
        feature_mask=np.zeros_like(wave, dtype=bool),
        min_pix_per_param=5,
    )

    # 8 pixels are insufficient for order 2 (15 required) or order 1 (10),
    # so the function falls back to the positive constant median continuum.
    np.testing.assert_allclose(fitted, np.nanmedian(flux))


def test_shape_mismatch_is_rejected():
    with pytest.raises(ValueError, match="identical shapes"):
        normSpecContinuum(
            "bad-shape",
            np.arange(10.0),
            np.ones(9),
            np.ones(10),
        )


def test_masked_absorption_equivalent_width_is_preserved():
    wave = np.linspace(5000.0, 5400.0, 4001)
    x = (wave - 5200.0) / 200.0
    continuum = 2.0 + 0.30 * x + 0.08 * x**2
    profile = 0.25 * np.exp(-0.5 * ((wave - 5200.0) / 3.0) ** 2)
    flux = continuum * (1.0 - profile)
    error = np.full_like(wave, 0.01)
    feature_mask = (wave >= 5185.0) & (wave <= 5215.0)

    flux_norm, _, fitted = normSpecContinuum(
        "ew-test",
        wave,
        flux,
        error,
        order=2,
        feature_mask=feature_mask,
    )

    ew_before = np.trapezoid(1.0 - flux / continuum, wave)
    ew_after = np.trapezoid(1.0 - flux_norm, wave)

    np.testing.assert_allclose(fitted, continuum, rtol=1e-10, atol=1e-10)
    assert ew_after == pytest.approx(ew_before, rel=1e-9, abs=1e-9)
