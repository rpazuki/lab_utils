import pandas as pd
import cobra

# Helper functions
def load_model(path):
    """Load a COBRA model from SBML or JSON."""
    if path.endswith(".xml") or path.endswith(".sbml"):
        return cobra.io.read_sbml_model(path)
    elif path.endswith(".json"):
        return cobra.io.load_json_model(path)
    else:
        raise ValueError("Unsupported model format")

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
              ) :
    df = df.copy()
    df[solution_column_name] = 0.0
    df[solution_status_column_name] = ""
    fluxes = []
    scale_factors = []
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
        else:
            fluxes.append([])

        scale_factors.append(scale_factor)
        #===============================
    fluxes_df = pd.concat([flux for flux in fluxes], axis=1, ignore_index=True)
    fluxes_df = fluxes_df.transpose()
    df["scaling_factor"] = scale_factors
    return df, fluxes_df

