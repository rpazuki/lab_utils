from pathlib import Path

import pandas as pd

from labUtils.utils import collate_by_strain, collate_by_strain_single_file


def _write_growth_csv(folder: Path, rows: list[dict]) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(folder / "growth_rates.csv", index=False)


def test_collate_by_strain_from_folders(tmp_path: Path):
    exp1 = tmp_path / "exp1"
    exp2 = tmp_path / "exp2"

    _write_growth_csv(
        exp1,
        [
            {"strain": "strainA1", "well": "A1", "value": 1.0},
            {"strain": "strainB1", "well": "B1", "value": 2.0},
        ],
    )
    _write_growth_csv(
        exp2,
        [
            {"strain": "strainA2", "well": "A2", "value": 3.0},
            {"strain": "strainB2", "well": "B2", "value": 4.0},
        ],
    )

    report = collate_by_strain(
        folders_list=[exp1, exp2],
        groupby_pattern=r"[A-Za-z]+",
        csv_output_file_name="report.csv",
    )

    assert len(report) == 2
    output_files = [Path(p) for p in report["output_file"]]
    assert all(path.exists() for path in output_files)

    output_by_group = {path.parent.name: pd.read_csv(path) for path in output_files}
    assert set(output_by_group) == {"strainA", "strainB"}
    assert set(output_by_group["strainA"]["strain"]) == {"strainA1", "strainA2"}
    assert set(output_by_group["strainA"]["experiment"]) == {"exp1", "exp2"}
    assert set(output_by_group["strainB"]["strain"]) == {"strainB1", "strainB2"}


def test_collate_by_strain_single_file(tmp_path: Path):
    input_csv = tmp_path / "single.csv"
    pd.DataFrame(
        [
            {"strain": "strainA1", "well": "A1", "value": 1.0},
            {"strain": "strainA2", "well": "A2", "value": 2.0},
            {"strain": "strainB1", "well": "B1", "value": 3.0},
        ]
    ).to_csv(input_csv, index=False)

    report = collate_by_strain_single_file(
        input_csv_file=input_csv,
        create_output_folder=False,
        groupby_pattern=r"[A-Za-z]+",
        csv_output_file_name="report.csv",
    )

    assert len(report) == 2
    output_files = [Path(p) for p in report["output_file"]]
    assert sorted(path.name for path in output_files) == ["strainA_report.csv", "strainB_report.csv"]

    strain_a_df = pd.read_csv(tmp_path / "strainA_report.csv")
    strain_b_df = pd.read_csv(tmp_path / "strainB_report.csv")
    assert set(strain_a_df["strain"]) == {"strainA1", "strainA2"}
    assert set(strain_b_df["strain"]) == {"strainB1"}
