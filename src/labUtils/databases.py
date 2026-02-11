
import pandas as pd
from pathlib import Path


def load_experiment_data(per_strain: bool = True,
                         replication: str = "replicates",
                         strain: str = "purB",
                         experiment: str = "mediabotJLF1",
                         well_column: str = "wells",
                         gr_column: str = "mv_mu_max",
                         od_cv_mean_threshold: float = 0.0,
                         od_cv_max_threshold: float = 0.0,
                         od_std_max_threshold: float = 0.0,
                         datasource_path:str = "H:/ROBOT_SCIENTIST/E_coli/Growth_rates/2025-10-31-27/processed",
                         levels_csv_file:str = "df_AMN_actual_medium_level.csv"
                         ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load experimental data, growth rates, and statistics; filter based on OD variability.
    """
    flg_has_statistic = True
    if per_strain:
        exp_data_path = Path(datasource_path) / f"{replication}_STRAINS" / strain / "AMN_dataset"
        growth_data_path = Path(datasource_path) / f"{replication}_STRAINS" / strain
        std_data_path = Path(datasource_path) / replication
        directories = [d for d in Path(std_data_path).iterdir()]
        std_data_list = [ pd.read_csv(d / "predictions.csv") for d in directories if (d / "predictions.csv").exists() ]
        for df , d in zip(std_data_list, directories):
            if "od600_mean" not in df.columns:
                flg_has_statistic = False
                break
            df["experiment"] = d.name
            # Coefficient of Variation (CV), also known as Relative Standard Deviation (RSD).
            df["od_CV"] = df['od600_std'] / df['od600_mean']

        if flg_has_statistic:
            std_data = pd.concat(std_data_list, ignore_index=True)
            std_data = std_data.groupby([well_column,"experiment"]).agg(
                od_mean_max = ( 'od600_mean', 'max' ),
                od_mean_mean = ( 'od600_mean', 'mean' ),
                od_std_max = ( 'od600_std', 'max' ),
                od_std_mean = ( 'od600_std', 'mean' ),
                od_cv_max = ('od_CV', 'max'),
                od_cv_mean = ('od_CV', 'mean'),
            ).reset_index()
        merging_directories = [d for d in (Path(datasource_path) / replication).iterdir() if not d.name.startswith(".")]
        df_parsed_data_list = [ pd.read_csv(d / "parsed_data.csv") for d in merging_directories if (d / "parsed_data.csv").exists() ]
        for df, d in zip(df_parsed_data_list, merging_directories):
            df["experiment"] = d.name
        df_parsed_data = pd.concat(df_parsed_data_list, ignore_index=True)
        df_parsed_data = df_parsed_data[df_parsed_data["strain"].str.startswith(strain)].copy().reset_index(drop=True)


    else:
        exp_data_path = Path(datasource_path) / replication / experiment / "AMN_dataset"
        growth_data_path = Path(datasource_path) / replication / experiment

        std_data = pd.read_csv(growth_data_path / "predictions.csv")
        if "od600_mean" not in std_data.columns:
            flg_has_statistic = False
            std_data = std_data.groupby(well_column).agg(
                od_mean_max = ( 'od600_mean', 'max' ),
                od_mean_mean = ( 'od600_mean', 'mean' ),
                od_std_max = ( 'od600_std', 'max' ),
                od_std_mean = ( 'od600_std', 'mean' ),
            ).reset_index()

            df_parsed_data = pd.read_csv(growth_data_path / "parsed_data.csv")

    #
    exp_data = pd.read_csv(exp_data_path / "df_flux.csv")
    growth_data = pd.read_csv(growth_data_path / "growth_rates.csv")
    df_levels = pd.read_csv(exp_data_path / levels_csv_file)

    # Filter growth_data to match the rows in exp_data
    # Since exp_data was created from growth_data where success=='ok', we need to apply the same filter
    # OR simply take the first len(exp_data) rows if they're already aligned
    growth_data = growth_data.reset_index(drop=True)
    growth_data = growth_data.loc[growth_data['success'], :]
    exp_data = exp_data.reset_index(drop=True)
    # Remove growth rate column from growth_data to avoid duplicates
    growth_data = growth_data.drop(columns=[gr_column])
    # Now join them - they should have the same number of rows
    combinded_data = pd.concat([exp_data.reset_index(drop=True), growth_data.reset_index(drop=True)],
                                axis=1)

    if flg_has_statistic:
        if per_strain:
            combinded_data = combinded_data.merge(std_data, on=["experiment", well_column], how="left")
        else:
            combinded_data = combinded_data.merge(std_data, on=[well_column], how="left")

    if flg_has_statistic and od_cv_mean_threshold > 0.0:
        combinded_data = combinded_data.loc[combinded_data["od_cv_mean"] <= od_cv_mean_threshold, :].copy().reset_index(drop=True)
        exp_data = combinded_data[exp_data.columns]
        if combinded_data.shape[0] < 2:
            print("Not enough data.")
            assert False
    elif flg_has_statistic and od_cv_max_threshold > 0.0:
        combinded_data = combinded_data.loc[combinded_data["od_cv_max"] <= od_cv_max_threshold, :].copy().reset_index(drop=True)
        exp_data = combinded_data[exp_data.columns]
        if combinded_data.shape[0] < 2:
            print("Not enough data.")
            assert False
    elif flg_has_statistic and od_std_max_threshold > 0.0:
        combinded_data = combinded_data.loc[combinded_data["od_std_max"] <= od_std_max_threshold, :].copy().reset_index(drop=True)
        exp_data = combinded_data[exp_data.columns]
        if combinded_data.shape[0] < 2:
            print("Not enough data.")
            assert False


    def update_uracil(val):
        for part in ["_0.5", "_22.4", "_200", "_2", "_8", "_64", "_640", "0"]:
                val = val.replace(part, "")
        return val
    combinded_data["supplements_unified"] = combinded_data["supplements"].apply(lambda x: update_uracil(str(x)) if pd.notna(x) else x)

    def classify(val):
        """Classify supplements into Sugar, Nucleo, Amino based on content"""
        if pd.isna(val) or str(val).lower() == 'nan' or str(val).strip() == '':
            return 'None'

        # Split by semicolon and convert to lowercase for comparison
        parts = str(val).lower().split(';')
        parts = [p.strip() for p in parts if p.strip()]

        if not parts:
            print(val)
            return 'None'

        # Define categories (all lowercase)
        sugars = ['glucose', 'succinate', 'sucrose', 'galactose', 'fructose', 'mannose',
                    'maltose', 'lactose', 'xylose', 'arabinose', 'ribose']

        nucleobases = ['adenine', 'uracil', 'guanine', 'cytosine', 'thymine']

        # Amino acids are typically 3 letters (all lowercase)
        amino_acids_3letter = ['ala', 'arg', 'asn', 'asp', 'cys', 'gln', 'glu', 'gly',
                                'his', 'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser',
                                'thr', 'trp', 'tyr', 'val']

        categories = []
        for part in parts:
            # Check if it's a sugar
            if any(sugar in part for sugar in sugars):
                if 'Sugar' not in categories:
                    categories.append('Sugar')
            # Check if it's a nucleobase
            elif any(nucleo in part for nucleo in nucleobases):
                if 'Nucleo' not in categories:
                    categories.append('Nucleo')
            # Check if it's a 3-letter amino acid or contains common amino acid patterns
            elif any(aa in part for aa in amino_acids_3letter) or len(part) == 3:
                if 'Amino' not in categories:
                    categories.append('Amino')

        if not categories:
            return 'Other'

        # Sort to ensure consistent ordering
        categories.sort()
        return '+'.join(categories)

    combinded_data["group"] = combinded_data["supplements_unified"].apply(classify)
    #===============================
    return combinded_data, exp_data, df_levels, df_parsed_data


def load_multiple_experiments_data(per_strain: bool = True,
                         replication: str = "replicates",
                         strains: list[str] = ["purB"],
                         experiments: list[str] = ["mediabotJLF1"],
                         well_column: str = "wells",
                         gr_column: str = "mv_mu_max",
                         od_cv_mean_threshold: float = 0.0,
                         od_cv_max_threshold: float = 0.0,
                         od_std_max_threshold: float = 0.0,
                         datasource_path:str = "H:/ROBOT_SCIENTIST/E_coli/Growth_rates/2025-10-31-27/processed",
                         levels_csv_file:str = "df_AMN_actual_medium_level.csv"
                         ) -> tuple[pd.DataFrame | None,
                                    pd.DataFrame | None,
                                    pd.DataFrame | None,
                                    pd.DataFrame | None]:
    """Load and combine data from multiple experiments."""
    df = None
    df_exp_data = None
    df_levels = None
    df_parsed_data = None
    if per_strain:
        for strain in strains:
            df_temp, df_exp_data_temp, df_levels_temp, df_pasred_data_temp = load_experiment_data(per_strain=per_strain,
                                            replication=replication,
                                            strain=strain,
                                            well_column=well_column,
                                            gr_column=gr_column,
                                            levels_csv_file="df_AMN_actual_medium_level.csv")
            if df is None:
                df = df_temp
                df_exp_data = df_exp_data_temp
                df_levels = df_levels_temp
                df_parsed_data = df_pasred_data_temp
            else:
                df = pd.concat([df, df_temp], ignore_index=True)
                df_exp_data = pd.concat([df_exp_data, df_exp_data_temp], ignore_index=True)
                df_parsed_data = pd.concat([df_parsed_data, df_pasred_data_temp], ignore_index=True)
    else:
        for experiment in experiments:
            df_temp, df_exp_data_temp, df_levels_temp, df_pasred_data_temp = load_experiment_data(per_strain=per_strain,
                                            replication=replication,
                                            experiment=experiment,
                                            well_column=well_column,
                                            gr_column=gr_column,
                                            levels_csv_file="df_AMN_actual_medium_level.csv")
            if df is None:
                df = df_temp
                df_exp_data = df_exp_data_temp
                df_levels = df_levels_temp
                df_parsed_data = df_pasred_data_temp
            else:
                df = pd.concat([df, df_temp], ignore_index=True)
                df_exp_data = pd.concat([df_exp_data, df_exp_data_temp], ignore_index=True)
                df_parsed_data = pd.concat([df_parsed_data, df_pasred_data_temp], ignore_index=True)
    return df, df_exp_data, df_levels, df_parsed_data


