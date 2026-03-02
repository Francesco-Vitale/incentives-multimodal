import amplpy
from amplpy import AMPL
import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path

for counter in range(6):

    def extract_entity(ampl, name, entity_type="variable"):
        """
        Extract an AMPL variable or parameter into a clean pandas DataFrame
        with a standardized value column equal to 'name'.
        """
        if entity_type == "variable":
            df = ampl.getVariable(name).getValues().toPandas()
        elif entity_type == "parameter":
            df = ampl.getParameter(name).getValues().toPandas()
        else:
            raise ValueError("entity_type must be 'variable' or 'parameter'")

        # Drop stray index column if present
        if 'index' in df.columns:
            df = df.drop(columns=['index'])

        # Identify value column (AMPL usually returns 1 value column)
        if len(df.columns) > 0:
            value_col = df.columns[-1]   # last column is safest choice
            df = df.rename(columns={value_col: name})
        else:
            df[name] = None

        return df

    BASE_DIR = Path(__file__).resolve().parent

    # Input directory
    input_dir = BASE_DIR / "input"

    if not os.path.exists(input_dir):
        print(f"Error: directory '{input_dir}' does not exist.")
        sys.exit(1)

    # Output directory
    output_dir = BASE_DIR / "./output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Network name
    network_name = (f"jatkasaari_network_{counter}")

    # Budget values for sensitivity analysis
    budgets = range(0,2_000_001,10_000)
    #budgets = range(0,50_001,5_000)
    budgets = (
	list(range(0, 11, 1)) +
        list(range(20, 101, 10)) +
        list(range(200, 1001, 100)) +
        list(range(2000, 10001, 1000)) +
        list(range(20000, 100001, 10000)) +
        list(range(200000, 1000001, 100000))
    )

    # Master DataFrame to collect results across budgets
    sensitivity_results = []

    # Create AMPL instance
    ampl = AMPL()

    try:
        # Read the model file (once, outside loop)
        model_path = os.path.join(input_dir, 'multim.mod')
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Model file not found: {model_path}")
        ampl.read(model_path)
        
        # Read the data file (once, outside loop; budgets will be overridden)
        data_path = os.path.join(input_dir, network_name + ".dat")
        if not os.path.exists(data_path):
            raise FileNotFoundError(f"Data file not found: {data_path}")
        ampl.read_data(data_path)
        
        # Solver setup: Gurobi for non-convex NLP/MPEC
        # Notes on Gurobi:
        # - This model is non-convex (exp, ^L power, complementarity), so use NonConvex=2 for global search.
        # - Complementarity (cond2) is treated as nonlinear; Gurobi handles MIQCP but may need spatial branch-and-bound.
        # - Settings: Time limit per solve (300s), threads (4), verbosity, and multi-start-like via pool.
        # - If convergence issues: Consider reformulating complementarity as big-M or SOS (e.g., via indicator constraints).
        solver_name = "gurobi"
        ampl.setOption('solver', solver_name)
        ampl.setOption("solver_msg", 0)
        ampl.setOption('gurobi_options', 
                    'NonConvex=2'  # Enable global optimization for non-convex problems
                    )  # Generate multiple solutions

        #print(f"Setup complete. Running sensitivity analysis for budgets: {budgets}")
        
        for i, B_val in enumerate(budgets):
            #print(f"\n--- Running solve {i+1}/{len(budgets)} with B = {B_val} ---")
            
            # Override budget parameter
            ampl.get_parameter("B").set(B_val)
            
            # Solve
            #print(f"Solving with {solver_name}...")
            ampl.solve()

            # Check status
            status = ampl.getValue('solve_result')
            #print(f"Solve status: {status}")

            if status in ['optimal', 'solved', 'feasible']:
                obj_val = ampl.getObjective('Total_Travel_Time').value()

                all_var_names = [name for name, _ in ampl.getVariables()]
                
                all_data = {}

                for var_name in all_var_names:
                    all_data[var_name] = extract_entity(ampl, var_name, "variable")

                all_data["ttpt"] = extract_entity(ampl, "ttpt", "parameter")

                od_variables = {}
                link_variables = {}

                for name, df in all_data.items():
                    # Example rule: separate by index structure
                    if df.index.nlevels == 1:
                        if len(df) == 30:
                            link_variables[name] = df
                        else:
                            od_variables[name] = df

                from functools import reduce

                if od_variables:
                    merged_od = reduce(
                        lambda left, right: left.merge(
                            right, left_index=True, right_index=True, how='outer'
                        ),
                        od_variables.values()
                    )

                merged_od['budget'] = B_val
                merged_od['objective'] = obj_val

                merged_od = merged_od.reset_index()
                merged_od.to_csv(f"./output/od_variables_{counter}.csv", index=False)

                sensitivity_results.append(merged_od)
            
                if link_variables:
                    merged_links = reduce(
                        lambda left, right: left.merge(
                            right, left_index=True, right_index=True, how='outer'
                        ),
                        link_variables.values()
                    )

                    merged_links = merged_links.reset_index()
                    merged_links.to_csv(f"./output/link_variables_{counter}.csv", index=False)

            else:
                #print(f"{solver_name} failed to converge for B={B_val}. Skipping...")
                # Optional: Add NaN row for this budget (dummy for single OD)
                dummy_row = pd.DataFrame({'budget': [B_val], 'qpt': [np.nan], 'qc': [np.nan], 'Ipt': [np.nan], 'ttpt': [np.nan], 'objective': [np.nan]}, index=['dummy'])
                sensitivity_results.append(dummy_row)

        # Combine all results into a single DataFrame
        if sensitivity_results:

            # Concatenate everything
            all_results = pd.concat(sensitivity_results, ignore_index=False)

            # Reset index safely and ensure OD column exists
            all_results = all_results.reset_index()
            if "OD" not in all_results.columns:
                all_results = all_results.rename(columns={"index": "OD"})

            # Convert everything except OD to numeric where possible
            for col in all_results.columns:
                if col != "OD":
                    try:
                        all_results[col] = pd.to_numeric(all_results[col])
                    except Exception:
                        pass

            # Save everything (no column filtering!)
            csv_path = os.path.join(output_dir, network_name + ".csv")
            all_results.to_csv(csv_path, index=False)

            #print("\nColumns stored:")
            #print(list(all_results.columns))
            '''
            # ---- Optional: automatic numeric summary by budget ----
            if "budget" in all_results.columns:

                numeric_cols = all_results.select_dtypes(include="number").columns.tolist()

                # Remove budget itself from summary targets
                numeric_cols = [c for c in numeric_cols if c != "budget"]

                if numeric_cols:
                    print("\nSummary Table (Budget vs. numeric averages):")
                    summary = (
                        all_results
                        .groupby("budget")[numeric_cols]
                        .mean()
                        .round(2)
                    )
                    print(summary)
                else:
                    print("No numeric columns available for summary.")
            else:
                print("No 'budget' column found. Skipping summary.")
            '''
        else:
            print("No results collected. Check solver convergence.")

    except FileNotFoundError as e:
        print(f"File error: {e}")
        print(f"Ensure 'multim.mod' and '{network_name}.dat' are in the input directory.")
    except Exception as e:
        print(f"Unexpected error: {e}")
        import traceback
        traceback.print_exc()  # Full stack trace for debugging
        sys.exit(1)
    finally:
        # Clean up
        print(f"Simulation {counter} completed.")
        del ampl
