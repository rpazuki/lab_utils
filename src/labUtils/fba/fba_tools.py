import pandas as pd
import cobra

from labUtils.databases import load_experiment_data
from labUtils.growth_rates import fit_max_growth_rate_per_series, fit_modified_gompertz_per_series, transform_to_log_n_n0, predict_modified_gompertz_per_series
from labUtils.media_bot import calculate_replicate_statistics_by_custom
from labUtils.utils import smart_join_drop_right

# Helper functions
def load_model(path):
    """Load a COBRA model from SBML or JSON."""
    if path.endswith(".xml") or path.endswith(".sbml"):
        return cobra.io.read_sbml_model(path)
    elif path.endswith(".json"):
        return cobra.io.load_json_model(path)
    else:
        raise ValueError("Unsupported model format")

def load_fba_data(per_strain: bool = True,
                replication: str = "replicates",
                strain: str = "purB",
                experiment: str = "mediabotJLF1",
                well_column: str = "wells",
                gr_column: str = "mv_mu_max",
                 group_cols: list[str] = ["well"],
                od_cv_mean_threshold: float = 0.0,
                od_cv_max_threshold: float = 0.0,
                od_std_max_threshold: float = 0.0,
                datasource_path:str = "H:/ROBOT_SCIENTIST/E_coli/Growth_rates/2025-10-31-27/processed",
                levels_csv_file:str = "df_AMN_actual_medium_level.csv",
                OD_0_averaging_window:int=4,
                moving_window_size:int=5,
                smoothing_iterations:int=2,
                smooth_window_size:int=2
                ) -> tuple[pd.DataFrame, list[str], list[str], list[str], pd.DataFrame]:
    df_data, _,df_levels, df_parsed_data = load_experiment_data(per_strain=per_strain,
                                                                replication=replication,
                                                                strain=strain,
                                                                experiment=experiment,
                                                                well_column=well_column,
                                                                gr_column=gr_column,
                                                                od_cv_mean_threshold=od_cv_mean_threshold,
                                                                od_cv_max_threshold=od_cv_max_threshold,
                                                                od_std_max_threshold=od_std_max_threshold,
                                                                datasource_path=datasource_path,
                                                                levels_csv_file=levels_csv_file)
    fixed_columns = []
    medium_columns = []
    supplement_columns = []

    level_index = -1
    max_value_index = -1

    for col in df_levels.columns:
        if col == "name":
            level_index = df_levels[col].values.tolist().index('level')
            max_value_index = df_levels[col].values.tolist().index('max_value')
            continue

        if df_levels[col].values[level_index] == 1 and df_levels[col].values[max_value_index] == 20.0:
            fixed_columns.append(col[:-2])  # Remove the trailing '_i' from the column name
        elif df_levels[col].values[level_index] == 1 and df_levels[col].values[max_value_index] != 20.0:
            medium_columns.append(col[:-2])  # Remove the trailing '_i' from the column name
        else:
            supplement_columns.append(col[:-2])  # Remove the trailing '_i' from the column name
    #
    if replication == "replicates":
        df_parsed_data = calculate_replicate_statistics_by_custom(df_parsed_data,
                                                                      strain_pattern="[A-Za-z]+",
                                                                      custom_rules={
                                                                          "SLAB":{ "direction": "alphabetical",
                                                                          "sample_size": 3 # e.g A1, A2, A3 - B1, B2, B3
                                                                          },
                                                                          "purB":{ "direction": "alphabetical",
                                                                                   "sample_size": 3 # e.g A1, A2, A3 - B1, B2, B3
                                                                          },
                                                                          "ilvI":{"direction": "alphabetical",
                                                                                  "sample_size": 3 # e.g A1, A2, A3 - B1, B2, B3
                                                                          },
                                                                          "WT":{"direction": "alphabetical",
                                                                                "sample_size": 2 # e.g A1, A2 - B1, B2
                                                                          },
                                                                          "BLANK":{"direction": "alphabetical",
                                                                                   "sample_size": 1 # e.g A1 - B1
                                                                          }})
    df_transformed = transform_to_log_n_n0(df_parsed_data,
                                           OD_0_averaging_window=OD_0_averaging_window,
                                           transformed_col="log_od_od0")
    df_fit_modified_gompertz = fit_modified_gompertz_per_series(df_transformed, value_col="log_od_od0")
    df_fit_max_growth_rate = fit_max_growth_rate_per_series(df_transformed, value_col="log_od_od0",
                                                            moving_window_size=moving_window_size,
                                                            smoothing_iterations=smoothing_iterations,
                                                            smooth_window_size=smooth_window_size)
    df_combined_fit = smart_join_drop_right(df_fit_modified_gompertz, df_fit_max_growth_rate)
    df_predict_modified_gompertz = predict_modified_gompertz_per_series(df_transformed, df_combined_fit)

    fit_cols = group_cols + ["y0", "A", "mu_max", "mv_mu_max", "lambda", "r2", "rmse", "n", "success", "message"]
    params_to_merge = df_combined_fit[[col for col in fit_cols if col in df_combined_fit.columns]]
    df_combined = df_predict_modified_gompertz.merge(params_to_merge, on=group_cols, how="left")



    return df_data, fixed_columns, medium_columns, supplement_columns, df_combined


def solve_fba(df:pd.DataFrame,
                     model:cobra.Model,
                     df_medium_columns:list,
                     scaled_columns:list[str] = [],
                     scale_factor:float = 1.0,
                     zero_uptake_columns:list[str] = [],
                     use_sfba:bool = False,
                     solution_column_name:str = "fba_growth",
                     solution_status_column_name:str = "fba_status",
              ) -> pd.DataFrame:
    df = df.copy()
    df[solution_column_name] = 0.0
    df[solution_status_column_name] = ""
    for _, row in df.iterrows():
        with model:
            # Set the medium according to the experimental data
            medium = model.medium.copy()
            for rxn_id in medium:
                medium[rxn_id] = 0.0  # Set high uptake rates for all exchange reactions
            for rxn_id in df_medium_columns:
                if rxn_id in scaled_columns:
                    medium[rxn_id] = row.get(rxn_id, 0.0) * scale_factor
                else:
                    medium[rxn_id] = row.get(rxn_id, 0.0)
            for rxn_id in zero_uptake_columns:
                medium[rxn_id] = 0.0
            model.medium = medium
            # Perform FBA
            if use_sfba:
                solution = cobra.flux_analysis.pfba(model)
            else:
                # model.objective = 'BIOMASS_Ec_iML1515_WT_75p37M'
                # model.objective = 'BIOMASS_Ec_iML1515_core_75p37M'
                solution = model.optimize()
            if solution.status == 'optimal':
                df.loc[df.index == row.name, solution_column_name] = solution.objective_value
                df.loc[df.index == row.name, solution_status_column_name] = solution.status
            else:
                df.loc[df.index == row.name, solution_status_column_name] = solution.status
        #===============================
    return df

def solve_fba_and_flux(df:pd.DataFrame,
                     model:cobra.Model,
                     df_medium_columns:list,
                     scaled_columns:list[str] = [],
                     scale_factor:float = 1.0,
                     zero_uptake_columns:list[str] = [],
                     use_sfba:bool = False,
                     solution_column_name:str = "fba_growth",
                     solution_status_column_name:str = "fba_status",
              ) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = df.copy()
    df[solution_column_name] = 0.0
    df[solution_status_column_name] = ""
    fluxes = []
    for _, row in df.iterrows():
        with model:
            # Set the medium according to the experimental data
            medium = model.medium.copy()
            for rxn_id in medium:
                medium[rxn_id] = 0.0  # Set high uptake rates for all exchange reactions
            for rxn_id in df_medium_columns:
                if rxn_id in scaled_columns:
                    medium[rxn_id] = row.get(rxn_id, 0.0) * scale_factor
                else:
                    medium[rxn_id] = row.get(rxn_id, 0.0)
            for rxn_id in zero_uptake_columns:
                medium[rxn_id] = 0.0
            model.medium = medium
            # Perform FBA
            if use_sfba:
                solution = cobra.flux_analysis.pfba(model)
            else:
                # model.objective = 'BIOMASS_Ec_iML1515_WT_75p37M'
                # model.objective = 'BIOMASS_Ec_iML1515_core_75p37M'
                solution = model.optimize()
            if solution.status == 'optimal':
                df.loc[df.index == row.name, solution_column_name] = solution.objective_value
                df.loc[df.index == row.name, solution_status_column_name] = solution.status
                fluxes.append(solution.fluxes.copy())
            else:
                df.loc[df.index == row.name, solution_status_column_name] = solution.status
                fluxes.append([])
        #===============================
    fluxes_df = pd.concat([flux for flux in fluxes], axis=1, ignore_index=True)
    fluxes_df = fluxes_df.transpose()
    return df, fluxes_df

def find_fba_solutions_scaling_factor(df:pd.DataFrame,
                     model:cobra.Model,
                     df_medium_columns:list,
                     scale_factor:float = 1.0,
                     zero_uptake_columns:list[str] = [],
                     use_sfba:bool = False,
                     max_iter:int = 20,
                     tol:float = 1e-6,
                     gr_column:str = "mv_mu_max",
                     solution_column_name:str = "fba_growth",
                     solution_status_column_name:str = "fba_status",
              ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = df.copy()
    df[solution_column_name] = 0.0
    df[solution_status_column_name] = ""
    fluxes = []
    scale_factors = []
    shadow_prices = []
    for _, row in df.iterrows():
        diff = float('inf')
        current_iteration = 0
        while diff > tol and current_iteration < max_iter:
            gr_value = row[gr_column]
            if diff == float('inf'):
                scale_factor = 1.0
            else:
                if solution.objective_value == 0:
                    print("fba_growth is zero ...")
                    fluxes.append([])
                    break
                scale_factor = scale_factor * (gr_value / solution.objective_value)
            with model:
                # Set the medium according to the experimental data
                medium = model.medium.copy()
                for rxn_id in medium:
                    medium[rxn_id] = 0.0  # Set high uptake rates for all exchange reactions
                for rxn_id in df_medium_columns:
                    medium[rxn_id] = row.get(rxn_id, 0.0) * scale_factor
                for rxn_id in zero_uptake_columns:
                    medium[rxn_id] = 0.0
                model.medium = medium
                # Perform FBA
                if use_sfba:
                    solution = cobra.flux_analysis.pfba(model)
                else:
                    # model.objective = 'BIOMASS_Ec_iML1515_WT_75p37M'
                    # model.objective = 'BIOMASS_Ec_iML1515_core_75p37M'
                    solution = model.optimize()
                if solution.status == 'optimal':
                    df.loc[df.index == row.name, solution_column_name] = solution.objective_value
                    df.loc[df.index == row.name, solution_status_column_name] = solution.status
                else:
                    df.loc[df.index == row.name, solution_status_column_name] = solution.status
                diff = abs(row[gr_column] - solution.objective_value)
            current_iteration += 1
        if solution.status == 'optimal':
            fluxes.append(solution.fluxes.copy())
            shadow_prices.append(solution.shadow_prices.copy()) # type: ignore
        else:
            fluxes.append([])
            shadow_prices.append([])

        scale_factors.append(scale_factor)
        #===============================
    fluxes_df = pd.concat([flux for flux in fluxes], axis=1, ignore_index=True)
    fluxes_df = fluxes_df.transpose()
    shadow_prices_df = pd.concat([sp for sp in shadow_prices], axis=1, ignore_index=True)
    shadow_prices_df = shadow_prices_df.transpose()
    df["scaling_factor"] = scale_factors
    return df, fluxes_df, shadow_prices_df

