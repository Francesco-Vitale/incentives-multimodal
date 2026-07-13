# Public Transport Incentive Optimization

This project implements a multimodal transportation optimization model that investigates how a limited incentive budget can encourage travelers to shift from private car travel to either public transport (PT) or a Park-and-Ride (P&R) alternative. The optimization model is formulated in **AMPL** and solved using the **Gurobi** optimizer. Python scripts are used to execute the optimization and generate summary figures.

## Project Structure

```text
.
├── input/
│   ├── jatkasaari_network.dat      # Network and demand data
│   └── multim.mod                  # AMPL optimization model
├── output/
│   ├── od_variables.csv
│   ├── link_variables.csv
│   ├── jatkasaari_network.csv
│   ├── jatkasaari_network_PT_perc.png
│   ├── jatkasaari_network_total_cost_benefit.png
│   ├── jatkasaari_network_TTT.png
│   └── jatkasaari_network_TTT_reduction.png
├── master_gurobi.py                # Runs the optimization model
├── plot_figures.py                 # Generates result figures
└── requirements.txt
```

## Requirements

The project requires:

- Python 3.10 or later
- AMPL
- Gurobi Optimizer (with a valid license)

Install the required Python packages using

```bash
pip install -r requirements.txt
```

or

```bash
pip install amplpy pandas numpy matplotlib
```

> **Note:** `amplpy` requires a working AMPL installation with access to the Gurobi solver.

## Workflow

### 1. Run the optimization

```bash
python master_gurobi.py
```

This script:

- loads the AMPL model and network data,
- performs a sensitivity analysis over the predefined incentive budgets,
- solves each optimization problem using Gurobi,
- stores the optimization results in the `output/` directory.

The generated CSV files are:

- `od_variables.csv`
- `link_variables.csv`
- `jatkasaari_network.csv`

### 2. Generate figures

```bash
python plot_figures.py
```

The script reads the generated CSV files and produces the following figures in the `output/` directory:

- `jatkasaari_network_PT_perc.png` – Public transport usage as a function of the incentive budget.
- `jatkasaari_network_total_cost_benefit.png` – Marginal benefit-to-cost ratio of increasing the incentive budget.
- `jatkasaari_network_TTT.png` – Total Travel Time (TTT) versus budget.
- `jatkasaari_network_TTT_reduction.png` – Percentage reduction in Total Travel Time relative to the baseline.

## Input Files

The `input/` directory contains:

- `multim.mod` – AMPL model formulation.
- `jatkasaari_network.dat` – Network topology, travel demand, and model parameters.

## Output Files

The optimization produces CSV files containing origin-destination variables, link-level variables, and aggregated simulation results. These outputs are subsequently used to generate summary figures illustrating changes in travel demand, public transport usage, network performance, and the effectiveness of the incentive budget.

## Methodology

The model evaluates how a fixed incentive budget can be allocated to encourage travelers to switch from private car travel to public transport or Park-and-Ride alternatives. A sensitivity analysis is performed across multiple budget levels, allowing the impact of increasing incentives on mode choice and overall network performance to be assessed.

The primary performance metric is **Total Travel Time (TTT)**, while additional outputs quantify public transport usage and the marginal benefit obtained from increasing the available incentive budget.