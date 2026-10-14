import pytest

from spectraPyle.runtime.runtime_adapter import adapt_gui_flat_to_schema
from spectraPyle.schema.schema import NormalizationConfig


def test_old_configs_get_continuum_defaults():
    cfg = NormalizationConfig.model_validate({"spectra_normalization": "median"})
    assert cfg.continuum_normalize is False
    assert cfg.continuum_norm_order == 2
    assert cfg.continuum_mask_features is True
    assert cfg.continuum_emission_mask_half_width == 12.0
    assert cfg.continuum_sigma_clip == 3.0


def test_flat_gui_adapter_does_not_inject_none_for_new_fields():
    structured = adapt_gui_flat_to_schema({"spectra_normalization": "median"})
    norm = structured["norm"]
    assert norm["continuum_normalize"] is False
    assert norm["continuum_norm_order"] == 2
    assert norm["continuum_mask_features"] is True
    assert norm["continuum_emission_mask_half_width"] == 12.0
    assert norm["continuum_sigma_clip"] == 3.0


def test_per_grism_polynomial_orders_are_supported():
    cfg = NormalizationConfig(
        spectra_normalization="median",
        continuum_normalize=True,
        continuum_norm_order={"arm1": 1, "arm2": 2, "default": 1},
    )
    assert cfg.continuum_norm_order["arm2"] == 2


def test_invalid_polynomial_order_is_rejected():
    with pytest.raises(ValueError, match="between 0 and 8"):
        NormalizationConfig(continuum_norm_order=9)


def test_template_and_continuum_normalization_are_incompatible():
    with pytest.raises(ValueError, match="cannot be combined"):
        NormalizationConfig(
            spectra_normalization="template",
            continuum_normalize=True,
        )


def test_continuum_normalized_output_is_dimensionless():
    pytest.importorskip("astropy")
    from spectraPyle.io.IO import output_units

    unit, scale = output_units(
        {
            "continuum_normalize": True,
            "spectra_normalization": "no_normalization",
            "conservation": "flux",
            "units": "erg/s/cm2/A",
            "units_scale_factor": 1e-16,
        }
    )
    assert unit == "dimensionless"
    assert scale == 1


def test_generic_yaml_template_validates():
    from pathlib import Path
    import yaml

    from spectraPyle.schema.schema import StackingConfig

    root = Path(__file__).resolve().parents[1]
    raw = yaml.safe_load((root / "configs/YAML/generic_template.yaml").read_text())
    cfg = StackingConfig.model_validate(raw)
    assert cfg.instrument.instrument_name == "generic"
    assert cfg.norm.lambda_norm_min == 5100.0
    assert cfg.norm.lambda_norm_max == 5500.0


def test_filename_records_continuum_normalization():
    from spectraPyle.io.filename_builder import build_filename

    base = {
        "survey_name": "generic",
        "grisms": ["arm1"],
        "data_release": "v1",
        "filename_in": "catalog.fits",
        "redshift_column_name": "z",
        "spectra_normalization": "median",
        "lambda_norm_rest": [None, None],
        "sigma_clipping_conditions": 3.0,
        "pixel_resampling": 2.0,
        "pixel_resampling_type": "lambda",
        "nyquist_sampling": 5.0,
        "continuum_mask_features": True,
        "continuum_emission_mask_half_width": 12.0,
        "continuum_sigma_clip": 3.0,
    }

    without = build_filename({**base, "continuum_normalize": False})
    with_cont = build_filename(
        {
            **base,
            "continuum_normalize": True,
            "continuum_norm_order": 2,
        }
    )

    assert "_cont2_" not in without
    assert "_cont2_" in with_cont
    assert without != with_cont
