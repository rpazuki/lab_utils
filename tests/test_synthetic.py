"""Tests for the in-silico synthetic dataset generator."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from labUtils.amn_mappings import MediumSource
from labUtils.synthetic import (
    EnumerationMode,
    EnumerationConfig,
    KineticsDefaultsConfig,
    PhenotypeMode,
    PhenotypeConfig,
    SupplementSpec,
    SyntheticDataset,
    attach_phenotype,
    build_flux_dataframe,
    enumerate_conditions,
    generate_community_dataset,
    generate_dataset,
    generate_dataset_from_args,
    generate_organism_dataset,
    write_dataset,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


GLUCOSE_MW = 180.156  # g/mol
URACIL_MW = 112.09
GLYCINE_MW = 75.07


def _mw_to_mmol(mass: float, mw: float) -> float:
    return (mass / mw) * 1000.0


@pytest.fixture
def simple_mappings_df() -> pd.DataFrame:
    """A minimal mappings_df bypassing `build_mappings` (so no PubChem network call)."""
    rows = [
        {
            "name": "glucose",
            "iupac_name": "",
            "other_names": "",
            "exchange_reaction": "EX_glc__D_e",
            "mass_per_litre": 4.0,
            "mmol_concentration": _mw_to_mmol(4.0, GLUCOSE_MW),
            "flux_upper_bound": 0.0,
            "source": MediumSource.SUPPLEMENT.value,
            "record_origin": "updated",
        },
        {
            "name": "uracil",
            "iupac_name": "",
            "other_names": "",
            "exchange_reaction": "EX_ura_e",
            "mass_per_litre": 0.0224,
            "mmol_concentration": _mw_to_mmol(0.0224, URACIL_MW),
            "flux_upper_bound": 0.0,
            "source": MediumSource.SUPPLEMENT.value,
            "record_origin": "updated",
        },
        {
            "name": "oxygen",
            "iupac_name": "",
            "other_names": "",
            "exchange_reaction": "EX_o2_e",
            "mass_per_litre": 0.0,
            "mmol_concentration": 0.0,
            "flux_upper_bound": 10.0,
            "source": MediumSource.FIXED.value,
            "record_origin": "updated",
        },
        {
            "name": "nh4",
            "iupac_name": "",
            "other_names": "",
            "exchange_reaction": "EX_nh4_e",
            "mass_per_litre": 1.0,
            "mmol_concentration": _mw_to_mmol(1.0, 18.04),
            "flux_upper_bound": 0.0,
            "source": MediumSource.MEDIUM.value,
            "record_origin": "updated",
        },
    ]
    return pd.DataFrame(rows)


@pytest.fixture
def toy_cobra_model():
    """Build a minimal in-memory cobra model with non-`EX_` and `EX_` exchange IDs."""
    cobra = pytest.importorskip("cobra")
    model = cobra.Model("toy")
    # Create metabolites
    glc_e = cobra.Metabolite("glc__D_e", compartment="e")
    biomass_c = cobra.Metabolite("biomass_c", compartment="c")
    o2_e = cobra.Metabolite("o2_e", compartment="e")

    # Reactions
    ex_glc = cobra.Reaction("EX_glc__D_e")
    ex_glc.lower_bound = -10
    ex_glc.upper_bound = 1000
    ex_glc.add_metabolites({glc_e: -1})

    ex_o2 = cobra.Reaction("EX_o2_e")
    ex_o2.lower_bound = -10
    ex_o2.upper_bound = 1000
    ex_o2.add_metabolites({o2_e: -1})

    biomass = cobra.Reaction("biomass")
    biomass.lower_bound = 0
    biomass.upper_bound = 1000
    biomass.add_metabolites({glc_e: -1, o2_e: -1, biomass_c: 1})

    sink = cobra.Reaction("DM_biomass")
    sink.lower_bound = 0
    sink.upper_bound = 1000
    sink.add_metabolites({biomass_c: -1})

    model.add_reactions([ex_glc, ex_o2, biomass, sink])
    model.objective = biomass
    return model


# ---------------------------------------------------------------------------
# Enumeration tests
# ---------------------------------------------------------------------------


class TestEnumeration:
    def test_cartesian_count(self):
        df = enumerate_conditions({
            "mode": "cartesian",
            "supplements": {
                "a": {"levels": [0.0, 1.0]},
                "b": {"levels": [0.0, 1.0, 2.0]},
                "c": {"levels": [0.0, 0.5, 1.0, 2.0]},
            },
        })
        assert len(df) == 2 * 3 * 4
        # Each combination unique
        assert df[["a", "b", "c"]].drop_duplicates().shape[0] == len(df)

    def test_custom_yaml_lists(self):
        df = enumerate_conditions({
            "mode": "custom",
            "supplements": {
                "a": {"levels": [0.1, 0.2]},
                "b": {"levels": [10, 20, 30]},
            },
        })
        assert len(df) == 2 * 3
        assert set(df["a"].unique()) == {0.1, 0.2}
        assert set(df["b"].unique()) == {10, 20, 30}

    def test_presence_absence_full_power_set(self):
        df = enumerate_conditions({
            "mode": "presence_absence",
            "supplements": {n: {"on_value": 1.0} for n in ["a", "b", "c", "d"]},
        })
        # 2^4 = 16
        assert len(df) == 16
        # All patterns unique
        bits = df[["a", "b", "c", "d"]].astype(int).drop_duplicates()
        assert bits.shape[0] == 16

    def test_presence_absence_bounded_k(self):
        df = enumerate_conditions({
            "mode": "presence_absence",
            "max_active": 2,
            "supplements": {n: {"on_value": 1.0} for n in ["a", "b", "c", "d"]},
        })
        # C(4,0) + C(4,1) + C(4,2) = 1 + 4 + 6 = 11
        assert len(df) == 11

    def test_random_lhs_seed_determinism(self):
        spec = {
            "mode": "random",
            "n_samples": 20,
            "seed": 123,
            "supplements": {
                "a": {"range": {"min": 0.0, "max": 1.0}},
                "b": {"range": {"min": -1.0, "max": 1.0}},
            },
        }
        df1 = enumerate_conditions(spec)
        df2 = enumerate_conditions(spec)
        pd.testing.assert_frame_equal(df1, df2)
        spec_diff = dict(spec, seed=999)
        df3 = enumerate_conditions(spec_diff)
        assert not df1["a"].equals(df3["a"])

    def test_random_lhs_n_samples_and_bounds(self):
        df = enumerate_conditions({
            "mode": "random",
            "n_samples": 50,
            "seed": 0,
            "supplements": {
                "a": {"range": {"min": 0.0, "max": 4.0}},
                "b": {"range": {"min": 1.0, "max": 2.0}},
            },
        })
        assert len(df) == 50
        assert df["a"].between(0.0, 4.0).all()
        assert df["b"].between(1.0, 2.0).all()

    def test_random_falls_back_when_scipy_missing(self, monkeypatch):
        import builtins

        real_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name.startswith("scipy"):
                raise ImportError("forced")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        df = enumerate_conditions({
            "mode": "random",
            "n_samples": 12,
            "seed": 7,
            "supplements": {"a": {"range": {"min": 0.0, "max": 1.0}}},
        })
        assert len(df) == 12

    def test_empty_supplements_returns_one_row(self):
        df = enumerate_conditions({"mode": "cartesian", "supplements": {}})
        assert len(df) == 1
        assert "condition_id" in df.columns
        assert df["supplements"].iloc[0] == ""

    def test_supplements_column_lists_active_only(self):
        df = enumerate_conditions({
            "mode": "cartesian",
            "supplements": {"a": {"levels": [0.0, 1.0]}, "b": {"levels": [0.0, 1.0]}},
        })
        # Row with only a active
        row_a = df[(df["a"] > 0) & (df["b"] == 0)].iloc[0]
        assert row_a["supplements"] == "a"
        row_all_off = df[(df["a"] == 0) & (df["b"] == 0)].iloc[0]
        assert row_all_off["supplements"] == ""

    def test_unknown_mode_raises(self):
        with pytest.raises(ValueError, match="Unknown enumeration mode"):
            enumerate_conditions({"mode": "bogus", "supplements": {"a": {"levels": [1.0]}}})


# ---------------------------------------------------------------------------
# Flux dataframe tests
# ---------------------------------------------------------------------------


class TestFluxDataframe:
    def test_columns_include_all_valid_exchanges_plus_mu_max(self, simple_mappings_df):
        conds = enumerate_conditions({
            "mode": "cartesian",
            "supplements": {"glucose": {"levels": [0.0, 4.0]}, "uracil": {"levels": [0.0]}},
        })
        flux = build_flux_dataframe(conds, simple_mappings_df)
        assert "mu_max" in flux.columns
        # All four exchanges should be columns (sorted, no suffix)
        exch_cols = [c for c in flux.columns if c != "mu_max"]
        assert exch_cols == sorted(["EX_glc__D_e", "EX_ura_e", "EX_o2_e", "EX_nh4_e"])

    def test_flux_formula_consistency(self, simple_mappings_df):
        conds = enumerate_conditions({
            "mode": "cartesian",
            "supplements": {"glucose": {"levels": [4.0]}},
        })
        kin = {"max_time": 24.0, "mv_mu_max_value": 0.5}
        flux = build_flux_dataframe(conds, simple_mappings_df, kinetics_defaults=kin)
        expected = _mw_to_mmol(4.0, GLUCOSE_MW) / (24.0 * 0.5)
        assert flux["EX_glc__D_e"].iloc[0] == pytest.approx(expected, rel=1e-6)

    def test_zero_concentration_yields_zero_flux(self, simple_mappings_df):
        conds = enumerate_conditions({
            "mode": "cartesian",
            "supplements": {"glucose": {"levels": [0.0]}},
        })
        flux = build_flux_dataframe(conds, simple_mappings_df)
        assert flux["EX_glc__D_e"].iloc[0] == 0.0

    def test_fixed_source_uses_flux_upper_bound(self, simple_mappings_df):
        conds = enumerate_conditions({"mode": "cartesian", "supplements": {}})
        flux = build_flux_dataframe(conds, simple_mappings_df)
        # FIXED source -> uses flux_upper_bound directly
        assert flux["EX_o2_e"].iloc[0] == 10.0

    def test_medium_source_uses_mapping_mmol(self, simple_mappings_df):
        conds = enumerate_conditions({"mode": "cartesian", "supplements": {}})
        kin = {"max_time": 24.0, "mv_mu_max_value": 0.5}
        flux = build_flux_dataframe(conds, simple_mappings_df, kinetics_defaults=kin)
        expected = _mw_to_mmol(1.0, 18.04) / (24.0 * 0.5)
        assert flux["EX_nh4_e"].iloc[0] == pytest.approx(expected, rel=1e-6)

    def test_exchange_suffix_applied(self, simple_mappings_df):
        conds = enumerate_conditions({
            "mode": "cartesian",
            "supplements": {"glucose": {"levels": [0.0]}},
        })
        flux = build_flux_dataframe(conds, simple_mappings_df, exchange_suffix="_i")
        assert "EX_glc__D_e_i" in flux.columns
        assert "EX_glc__D_e" not in flux.columns
        assert "mu_max" in flux.columns
        assert "mu_max_i" not in flux.columns

    def test_no_EX_string_in_synthetic_module(self):
        """Regression: synthetic module must not hard-code the BiGG `EX_` prefix."""
        src = Path(__file__).resolve().parent.parent / "src" / "labUtils" / "synthetic.py"
        text = src.read_text(encoding="utf-8")
        assert "EX_" not in text, (
            "synthetic.py must not hard-code 'EX_' anywhere — exchange IDs come from "
            "the organism's mapping (which is grounded in its SBML)."
        )

    def test_exchange_ids_non_EX_prefix_handled(self):
        """Exchange IDs that do NOT start with 'EX_' must still flow through."""
        mappings = pd.DataFrame([
            {
                "name": "glucose", "iupac_name": "", "other_names": "",
                "exchange_reaction": "R_EX_glc_e",
                "mass_per_litre": 4.0,
                "mmol_concentration": _mw_to_mmol(4.0, GLUCOSE_MW),
                "flux_upper_bound": 0.0,
                "source": MediumSource.SUPPLEMENT.value,
                "record_origin": "updated",
            },
            {
                "name": "glycine", "iupac_name": "", "other_names": "",
                "exchange_reaction": "sink_glycine_c",  # KEGG-style sink, no EX_
                "mass_per_litre": 0.0, "mmol_concentration": 0.0,
                "flux_upper_bound": 5.0,
                "source": MediumSource.FIXED.value,
                "record_origin": "updated",
            },
        ])
        conds = enumerate_conditions({
            "mode": "cartesian",
            "supplements": {"glucose": {"levels": [4.0]}},
        })
        flux = build_flux_dataframe(conds, mappings)
        assert "R_EX_glc_e" in flux.columns
        assert "sink_glycine_c" in flux.columns
        assert flux["sink_glycine_c"].iloc[0] == 5.0

    def test_supplement_overrides_medium_precedence_per_row(self):
        """When a supplement and medium share an exchange, the supplement value wins per row."""
        mappings = pd.DataFrame([
            {
                "name": "glucose", "iupac_name": "", "other_names": "",
                "exchange_reaction": "EX_glc_e",
                "mass_per_litre": 4.0,
                "mmol_concentration": _mw_to_mmol(4.0, GLUCOSE_MW),
                "flux_upper_bound": 0.0,
                "source": MediumSource.SUPPLEMENT.value,
                "record_origin": "updated",
            },
            {
                "name": "glucose_in_medium", "iupac_name": "", "other_names": "",
                "exchange_reaction": "EX_glc_e",
                "mass_per_litre": 0.5,
                "mmol_concentration": _mw_to_mmol(0.5, GLUCOSE_MW),
                "flux_upper_bound": 0.0,
                "source": MediumSource.MEDIUM.value,
                "record_origin": "updated",
            },
        ])
        conds = enumerate_conditions({
            "mode": "cartesian",
            "supplements": {"glucose": {"levels": [0.0, 4.0]}},
        })
        flux = build_flux_dataframe(conds, mappings)
        # First row: glucose=0 -> falls back to medium value
        zero_row = conds[conds["glucose"] == 0.0].index[0]
        on_row = conds[conds["glucose"] == 4.0].index[0]
        # When glucose=0, supplement override does not kick in; medium value is set in template.
        assert flux.loc[zero_row, "EX_glc_e"] == pytest.approx(
            _mw_to_mmol(0.5, GLUCOSE_MW) / (24.0 * 0.5), rel=1e-6
        )
        # When glucose=4, supplement overrides.
        assert flux.loc[on_row, "EX_glc_e"] == pytest.approx(
            _mw_to_mmol(4.0, GLUCOSE_MW) / (24.0 * 0.5), rel=1e-6
        )

    def test_supplement_mmol_override_takes_precedence_over_mass(self, simple_mappings_df):
        conds = enumerate_conditions({
            "mode": "cartesian",
            "supplements": {"glucose": {"levels": [4.0]}},
        })
        flux = build_flux_dataframe(
            conds,
            simple_mappings_df,
            supplement_mmol_overrides={"glucose": 10.0},
        )
        expected = 10.0 / (24.0 * 0.5)
        assert flux["EX_glc__D_e"].iloc[0] == pytest.approx(expected, rel=1e-6)


# ---------------------------------------------------------------------------
# Phenotype attachment tests
# ---------------------------------------------------------------------------


class TestPhenotype:
    def test_empty_mode_leaves_mu_max_zero(self, simple_mappings_df):
        conds = enumerate_conditions({
            "mode": "cartesian", "supplements": {"glucose": {"levels": [0.0, 4.0]}},
        })
        flux = build_flux_dataframe(conds, simple_mappings_df)
        out = attach_phenotype(flux, conds, mode="empty")
        assert (out["mu_max"] == 0.0).all()

    def test_formula_evaluation(self, simple_mappings_df):
        conds = enumerate_conditions({
            "mode": "cartesian", "supplements": {"glucose": {"levels": [0.0, 4.0]}},
        })
        flux = build_flux_dataframe(conds, simple_mappings_df)
        out = attach_phenotype(flux, conds, mode="formula", formula="0.6 * (glucose > 0)")
        # Where glucose > 0: mu_max == 0.6; else 0.0
        for i, row in conds.iterrows():
            expected = 0.6 if row["glucose"] > 0 else 0.0
            assert out.loc[i, "mu_max"] == pytest.approx(expected)

    def test_formula_rejects_disallowed_names(self, simple_mappings_df):
        conds = enumerate_conditions({
            "mode": "cartesian", "supplements": {"glucose": {"levels": [0.0]}},
        })
        flux = build_flux_dataframe(conds, simple_mappings_df)
        with pytest.raises(ValueError, match="disallowed name"):
            attach_phenotype(flux, conds, mode="formula", formula="__import__('os')")

    def test_formula_rejects_attribute_access(self, simple_mappings_df):
        conds = enumerate_conditions({
            "mode": "cartesian", "supplements": {"glucose": {"levels": [0.0]}},
        })
        flux = build_flux_dataframe(conds, simple_mappings_df)
        with pytest.raises(ValueError, match="Attribute access"):
            attach_phenotype(flux, conds, mode="formula", formula="glucose.real")

    def test_formula_noise_deterministic(self, simple_mappings_df):
        conds = enumerate_conditions({
            "mode": "cartesian", "supplements": {"glucose": {"levels": [0.0, 4.0]}},
        })
        flux = build_flux_dataframe(conds, simple_mappings_df)
        a = attach_phenotype(flux, conds, mode="formula", formula="0.5", noise_std=0.1, seed=11)
        b = attach_phenotype(flux, conds, mode="formula", formula="0.5", noise_std=0.1, seed=11)
        pd.testing.assert_series_equal(a["mu_max"], b["mu_max"])

    def test_fba_mode_matches_solve_fba(self, toy_cobra_model):
        # Tiny mapping for the toy model
        mappings = pd.DataFrame([
            {
                "name": "glucose", "iupac_name": "", "other_names": "",
                "exchange_reaction": "EX_glc__D_e",
                "mass_per_litre": 4.0,
                "mmol_concentration": _mw_to_mmol(4.0, GLUCOSE_MW),
                "flux_upper_bound": 0.0,
                "source": MediumSource.SUPPLEMENT.value,
                "record_origin": "updated",
            },
            {
                "name": "oxygen", "iupac_name": "", "other_names": "",
                "exchange_reaction": "EX_o2_e",
                "mass_per_litre": 0.0, "mmol_concentration": 0.0,
                "flux_upper_bound": 10.0,
                "source": MediumSource.FIXED.value,
                "record_origin": "updated",
            },
        ])
        conds = enumerate_conditions({
            "mode": "cartesian",
            "supplements": {"glucose": {"levels": [4.0]}},
        })
        flux = build_flux_dataframe(conds, mappings)
        medium_cols = [c for c in flux.columns if c != "mu_max"]
        out = attach_phenotype(flux, conds, mode="fba", model=toy_cobra_model, medium_columns=medium_cols)
        # Toy model: max growth is min(glucose_uptake, o2_uptake) -> bounded by min flux
        assert out["mu_max"].iloc[0] > 0.0

    def test_fba_mode_requires_model_and_medium(self, simple_mappings_df):
        conds = enumerate_conditions({
            "mode": "cartesian", "supplements": {"glucose": {"levels": [0.0]}},
        })
        flux = build_flux_dataframe(conds, simple_mappings_df)
        with pytest.raises(ValueError, match="requires `model`"):
            attach_phenotype(flux, conds, mode="fba")

    def test_unknown_phenotype_mode_raises(self, simple_mappings_df):
        conds = enumerate_conditions({
            "mode": "cartesian", "supplements": {"glucose": {"levels": [0.0]}},
        })
        flux = build_flux_dataframe(conds, simple_mappings_df)
        with pytest.raises(ValueError, match="Unknown phenotype mode"):
            attach_phenotype(flux, conds, mode="bogus")


# ---------------------------------------------------------------------------
# generate_dataset / community / write
# ---------------------------------------------------------------------------


def _community_config_dict() -> dict:
    """Two organisms with different exchange-ID conventions for the same metabolite."""
    return {
        "synthetic": {
            "organisms": {
                "A": {
                    "supplement_to_exchange_map": {"glucose": "EX_glc__D_e"},
                    "molecular_weights": {"glucose": GLUCOSE_MW},
                },
                "B": {
                    "supplement_to_exchange_map": {"glucose": "R_EX_glc_e"},
                    "exchange_suffix": "_i",
                    "molecular_weights": {"glucose": GLUCOSE_MW},
                },
            },
            "enumeration": {
                "mode": "cartesian",
                "supplements": {"glucose": {"levels": [0.0, 4.0]}},
            },
            "phenotype": {"mode": "empty"},
            "growth_kinetics_defaults": {"max_time": 24.0, "mv_mu_max_value": 0.5},
        }
    }


class TestGenerateDataset:
    def test_single_organism_basic(self):
        cfg = {
            "synthetic": {
                "organisms": {"ecoli": {"supplement_to_exchange_map": {"glucose": "EX_glc__D_e"}}},
                "enumeration": {
                    "mode": "cartesian",
                    "supplements": {"glucose": {"levels": [0.0, 4.0]}},
                },
                "phenotype": {"mode": "empty"},
            }
        }
        ds = generate_dataset(cfg)
        assert isinstance(ds, SyntheticDataset)
        assert list(ds.flux_tables.keys()) == ["ecoli"]
        # flux_df shortcut works
        assert "mu_max" in ds.flux_df.columns
        # conditions_df has 2 rows
        assert len(ds.conditions_df) == 2

    def test_community_alignment(self):
        ds = generate_dataset(_community_config_dict())
        assert set(ds.flux_tables.keys()) == {"A", "B"}
        # Same number of rows across all tables
        n = len(ds.conditions_df)
        assert n == 2
        assert len(ds.medium_df) == n
        assert len(ds.flux_tables["A"]) == n
        assert len(ds.flux_tables["B"]) == n
        # B has suffix; A does not
        assert any(c.endswith("_i") for c in ds.flux_tables["B"].columns if c != "mu_max")
        assert all(not c.endswith("_i") for c in ds.flux_tables["A"].columns if c != "mu_max")

    def test_community_medium_is_keyed_by_common_name(self):
        ds = generate_dataset(_community_config_dict())
        # Common-name key 'glucose' should be in medium_df, not an exchange ID
        assert "glucose" in ds.medium_df.columns
        assert "EX_glc__D_e" not in ds.medium_df.columns
        assert "R_EX_glc_e" not in ds.medium_df.columns

    def test_community_exchange_id_mismatch_resolved(self):
        ds = generate_dataset(_community_config_dict())
        # Organism A uses EX_glc__D_e; organism B uses R_EX_glc_e
        a_cols = set(ds.flux_tables["A"].columns)
        b_cols = set(ds.flux_tables["B"].columns)
        assert "EX_glc__D_e" in a_cols
        assert "R_EX_glc_e_i" in b_cols
        # Both should have a non-zero glucose flux on the glucose=4.0 row
        on_idx = ds.conditions_df[ds.conditions_df["glucose"] > 0].index[0]
        assert ds.flux_tables["A"].loc[on_idx, "EX_glc__D_e"] > 0
        assert ds.flux_tables["B"].loc[on_idx, "R_EX_glc_e_i"] > 0

    def test_missing_organisms_raises(self):
        with pytest.raises(ValueError, match="at least one organism"):
            generate_dataset({"synthetic": {"enumeration": {"mode": "cartesian", "supplements": {}}}})

    def test_generate_dataset_from_args_matches_config_path(self):
        cfg = {
            "synthetic": {
                "organisms": {
                    "ecoli": {
                        "supplement_to_exchange_map": {"glucose": "EX_glc__D_e"},
                        "molecular_weights": {"glucose": GLUCOSE_MW},
                    }
                },
                "enumeration": {
                    "mode": "cartesian",
                    "supplements": {"glucose": {"levels": [0.0, 4.0]}},
                },
                "phenotype": {"mode": "empty"},
                "growth_kinetics_defaults": {"max_time": 24.0, "mv_mu_max_value": 0.5},
            }
        }

        from_cfg = generate_dataset(cfg)
        from_args = generate_dataset_from_args(
            organisms=cfg["synthetic"]["organisms"],
            enumeration=cfg["synthetic"]["enumeration"],
            phenotype=cfg["synthetic"]["phenotype"],
            kinetics_defaults=cfg["synthetic"]["growth_kinetics_defaults"],
        )

        pd.testing.assert_frame_equal(from_cfg.conditions_df, from_args.conditions_df)
        pd.testing.assert_frame_equal(from_cfg.medium_df, from_args.medium_df)
        pd.testing.assert_frame_equal(from_cfg.flux_df, from_args.flux_df)

    def test_generate_organism_dataset_matches_config_path(self):
        cfg = {
            "synthetic": {
                "organisms": {
                    "organism": {
                        "sbml": "organism.xml",
                        "supplement_to_exchange_map": {"glucose": "EX_glc__D_e"},
                        "exchange_suffix": "_i",
                    }
                },
                "enumeration": {
                    "mode": "cartesian",
                    "supplements": {"glucose": {"levels": [0.0, 4.0]}},
                },
                "phenotype": {"mode": "empty"},
            }
        }

        from_cfg = generate_dataset(cfg)
        from_args = generate_organism_dataset(
            organism_sbml_path="organism.xml",
            supplement_to_exchange_map={"glucose": "EX_glc__D_e"},
            exchange_suffix="_i",
            enumeration=EnumerationConfig(
                mode=EnumerationMode.CARTESIAN,
                supplements={"glucose": SupplementSpec(levels=[0.0, 4.0])},
            ),
            phenotype=PhenotypeConfig(mode=PhenotypeMode.EMPTY),
        )

        pd.testing.assert_frame_equal(from_cfg.conditions_df, from_args.conditions_df)
        pd.testing.assert_frame_equal(from_cfg.medium_df, from_args.medium_df)
        pd.testing.assert_frame_equal(from_cfg.flux_df, from_args.flux_df)

    def test_generate_organism_dataset_requires_enumeration(self):
        with pytest.raises(ValueError, match="`enumeration` is required"):
            generate_organism_dataset(
                supplement_to_exchange_map={"glucose": "EX_glc__D_e"},
                enumeration=None,
            )

    def test_generate_organism_dataset_requires_mapping(self):
        with pytest.raises(ValueError, match="non-empty mapping"):
            generate_organism_dataset(
                supplement_to_exchange_map={},
                enumeration=EnumerationConfig(
                    mode=EnumerationMode.CARTESIAN,
                    supplements={"glucose": SupplementSpec(levels=[0.0, 4.0])},
                ),
            )

    def test_generate_community_dataset_matches_config_path(self):
        cfg = {
            "synthetic": {
                "organisms": {
                    "ecoli": {
                        "sbml": "ecoli.xml",
                        "supplement_to_exchange_map": {"glucose": "EX_glc__D_e"},
                        "exchange_suffix": "_i",
                    },
                    "bsubtilis": {
                        "sbml": "bsubtilis.xml",
                        "supplement_to_exchange_map": {"glucose": "EX_glc__D_e"},
                        "exchange_suffix": "_i",
                    },
                },
                "enumeration": {
                    "mode": "cartesian",
                    "supplements": {"glucose": {"levels": [0.0, 4.0]}},
                },
                "phenotype": {"mode": "empty"},
            }
        }

        from_cfg = generate_dataset(cfg)
        from_args = generate_community_dataset(
            community_sbml_paths={"ecoli": "ecoli.xml", "bsubtilis": "bsubtilis.xml"},
            supplement_to_exchange_map={"glucose": "EX_glc__D_e"},
            exchange_suffix="_i",
            enumeration=EnumerationConfig(
                mode=EnumerationMode.CARTESIAN,
                supplements={"glucose": SupplementSpec(levels=[0.0, 4.0])},
            ),
            phenotype=PhenotypeConfig(mode=PhenotypeMode.EMPTY),
        )

        pd.testing.assert_frame_equal(from_cfg.conditions_df, from_args.conditions_df)
        pd.testing.assert_frame_equal(from_cfg.medium_df, from_args.medium_df)
        assert set(from_args.flux_tables) == {"ecoli", "bsubtilis"}

    def test_generate_community_dataset_requires_sbml_paths(self):
        with pytest.raises(ValueError, match="non-empty mapping of organism name -> SBML path"):
            generate_community_dataset(
                community_sbml_paths={},
                supplement_to_exchange_map={"glucose": "EX_glc__D_e"},
                enumeration=EnumerationConfig(
                    mode=EnumerationMode.CARTESIAN,
                    supplements={"glucose": SupplementSpec(levels=[0.0, 4.0])},
                ),
            )

    def test_generate_organism_dataset_requires_sbml_for_fba(self):
        with pytest.raises(ValueError, match="requires an SBML path for organism 'organism'"):
            generate_organism_dataset(
                supplement_to_exchange_map={"glucose": "EX_glc__D_e"},
                enumeration=EnumerationConfig(
                    mode=EnumerationMode.CARTESIAN,
                    supplements={"glucose": SupplementSpec(levels=[0.0, 4.0])},
                ),
                phenotype=PhenotypeConfig(mode=PhenotypeMode.FBA),
            )

    def test_supplement_dataclass_roundtrip_keeps_dict_shape(self):
        cfg = EnumerationConfig(
            mode=EnumerationMode.CARTESIAN,
            supplements={
                "glucose": SupplementSpec(levels=[0.0, 2.0, 4.0]),
                "uracil": SupplementSpec(binary=True, on_value=0.0224),
                "glycine": SupplementSpec(range_min=0.0, range_max=0.5, range_n=5),
            },
        )
        as_dict = cfg.to_dict()
        assert isinstance(as_dict["supplements"], dict)
        assert set(as_dict["supplements"].keys()) == {"glucose", "uracil", "glycine"}
        assert as_dict["supplements"]["glucose"]["levels"] == [0.0, 2.0, 4.0]

    def test_enumeration_override_works_when_yaml_section_missing(self):
        cfg = {
            "synthetic": {
                "organisms": {
                    "ecoli": {
                        "supplement_to_exchange_map": {"glucose": "EX_glc__D_e"},
                        "molecular_weights": {"glucose": GLUCOSE_MW},
                    }
                }
            }
        }
        enum = EnumerationConfig(
            mode=EnumerationMode.CUSTOM,
            supplements={"glucose": SupplementSpec(levels=[1.5, 2.5])},
        )
        ds = generate_dataset(cfg, enumeration=enum)
        assert set(ds.conditions_df["glucose"].astype(float).tolist()) == {1.5, 2.5}

    def test_phenotype_override_precedence_over_yaml(self):
        cfg = {
            "synthetic": {
                "organisms": {
                    "ecoli": {
                        "supplement_to_exchange_map": {"glucose": "EX_glc__D_e"},
                        "molecular_weights": {"glucose": GLUCOSE_MW},
                    }
                },
                "enumeration": {
                    "mode": "cartesian",
                    "supplements": {"glucose": {"levels": [0.0, 4.0]}},
                },
                "phenotype": {"mode": "empty"},
            }
        }
        ds = generate_dataset(cfg, phenotype=PhenotypeConfig(mode=PhenotypeMode.FORMULA, formula="0.25"))
        assert np.allclose(ds.flux_df["mu_max"].to_numpy(dtype=float), 0.25)

    def test_kinetics_defaults_override_precedence_over_yaml(self):
        cfg = {
            "synthetic": {
                "organisms": {
                    "ecoli": {
                        "supplement_to_exchange_map": {"glucose": "EX_glc__D_e"},
                        "molecular_weights": {"glucose": GLUCOSE_MW},
                    }
                },
                "enumeration": {
                    "mode": "cartesian",
                    "supplements": {"glucose": {"levels": [4.0]}},
                },
                "growth_kinetics_defaults": {"max_time": 24.0, "mv_mu_max_value": 0.5},
            }
        }
        baseline = generate_dataset(cfg)
        overridden = generate_dataset(
            cfg,
            kinetics_defaults=KineticsDefaultsConfig(max_time=12.0, mv_mu_max_value=0.5),
        )
        base_flux = float(baseline.flux_df["EX_glc__D_e"].iloc[0])
        override_flux = float(overridden.flux_df["EX_glc__D_e"].iloc[0])
        assert override_flux == pytest.approx(base_flux * 2.0, rel=1e-6)

    def test_enumeration_mmol_per_liter_override_applies_to_supplement_flux(self):
        cfg = {
            "synthetic": {
                "organisms": {
                    "ecoli": {
                        "supplement_to_exchange_map": {"glucose": "EX_glc__D_e"},
                        "molecular_weights": {"glucose": GLUCOSE_MW},
                    }
                },
                "enumeration": {
                    "mode": "cartesian",
                    "supplements": {
                        "glucose": {
                            "levels": [4.0],
                            "mmol_per_liter": 5.0,
                        }
                    },
                },
            }
        }
        ds = generate_dataset(cfg)
        expected = 5.0 / (24.0 * 0.5)
        assert float(ds.flux_df["EX_glc__D_e"].iloc[0]) == pytest.approx(expected, rel=1e-6)

    def test_enumeration_mmol_per_liter_wins_over_mol_per_liter(self):
        cfg = {
            "synthetic": {
                "organisms": {
                    "ecoli": {
                        "supplement_to_exchange_map": {"glucose": "EX_glc__D_e"},
                        "molecular_weights": {"glucose": GLUCOSE_MW},
                    }
                },
                "enumeration": {
                    "mode": "cartesian",
                    "supplements": {
                        "glucose": {
                            "levels": [4.0],
                            "mmol_per_liter": 3.0,
                            "mol_per_liter": 0.2,
                        }
                    },
                },
            }
        }
        ds = generate_dataset(cfg)
        expected = 3.0 / (24.0 * 0.5)
        assert float(ds.flux_df["EX_glc__D_e"].iloc[0]) == pytest.approx(expected, rel=1e-6)


class TestWriteDataset:
    def test_write_and_reload_roundtrip(self, tmp_path: Path):
        cfg = {
            "synthetic": {
                "organisms": {"ecoli": {"supplement_to_exchange_map": {"glucose": "EX_glc__D_e"}}},
                "enumeration": {
                    "mode": "cartesian",
                    "supplements": {"glucose": {"levels": [0.0, 4.0]}},
                },
                "phenotype": {"mode": "empty"},
            }
        }
        ds = generate_dataset(cfg)
        paths = write_dataset(ds, tmp_path)
        assert paths["conditions"].exists()
        assert paths["medium"].exists()
        assert paths["ecoli"].exists()
        # Reload and compare
        flux_round = pd.read_csv(paths["ecoli"])
        pd.testing.assert_frame_equal(flux_round.reset_index(drop=True), ds.flux_df.reset_index(drop=True))


# ---------------------------------------------------------------------------
# End-to-end with toy COBRA model
# ---------------------------------------------------------------------------


def test_end_to_end_solve_fba_on_generated_flux(toy_cobra_model):
    from labUtils.fba.fba_tools import solve_fba

    mappings = pd.DataFrame([
        {
            "name": "glucose", "iupac_name": "", "other_names": "",
            "exchange_reaction": "EX_glc__D_e",
            "mass_per_litre": 4.0,
            "mmol_concentration": _mw_to_mmol(4.0, GLUCOSE_MW),
            "flux_upper_bound": 0.0,
            "source": MediumSource.SUPPLEMENT.value,
            "record_origin": "updated",
        },
        {
            "name": "oxygen", "iupac_name": "", "other_names": "",
            "exchange_reaction": "EX_o2_e",
            "mass_per_litre": 0.0, "mmol_concentration": 0.0,
            "flux_upper_bound": 10.0,
            "source": MediumSource.FIXED.value,
            "record_origin": "updated",
        },
    ])
    conds = enumerate_conditions({
        "mode": "cartesian",
        "supplements": {"glucose": {"levels": [0.0, 4.0]}},
    })
    flux = build_flux_dataframe(conds, mappings)
    medium_cols = [c for c in flux.columns if c != "mu_max"]
    solved = solve_fba(flux, toy_cobra_model, medium_cols)
    # Row with glucose > 0 should be optimal with non-zero growth
    on_idx = conds[conds["glucose"] > 0].index[0]
    assert solved.loc[on_idx, "fba_status"] == "optimal"
    assert solved.loc[on_idx, "fba_growth"] > 0
