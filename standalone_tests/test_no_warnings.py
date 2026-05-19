"""Tests for warning behavior in fuzzy supplement matching."""

import logging

import pandas as pd

from labUtils.amn_mappings import build_mappings, load_default_iml1515_mapping


def _fuzzy_growth_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "well": ["A1"],
            "supplements": ["TotallyUnknownSupplement"],
            "mu_max": [0.2],
        }
    )


def test_verbose_false_suppresses_fuzzy_warnings(caplog):
    caplog.set_level(logging.WARNING)

    build_mappings(
        _fuzzy_growth_df(),
        supplement_to_exchange_map=load_default_iml1515_mapping(),
        fuzzy_threshold=0.6,
        verbose=False,
    )

    assert not any("fuzzy match found" in rec.message.lower() for rec in caplog.records)


def test_verbose_true_emits_fuzzy_warnings(caplog):
    caplog.set_level(logging.WARNING)

    build_mappings(
        _fuzzy_growth_df(),
        supplement_to_exchange_map=load_default_iml1515_mapping(),
        fuzzy_threshold=0.6,
        verbose=True,
    )

    assert any("could not found a match" in rec.message.lower() for rec in caplog.records)
