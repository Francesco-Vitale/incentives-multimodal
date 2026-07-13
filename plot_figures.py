import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.ticker import ScalarFormatter

# -------------------------------------------------------------
# USER SETTINGS
# -------------------------------------------------------------

output_dir = r".\output"
saveFig = True

variables_to_plot = [
    "qpt", "qc",
    "qpt_share", "qc_share",
    "PT_TT", "Car_TT",
    "Upt", "Ucar",
    "objective"
]

# -------------------------------------------------------------


def main():

    csv_path = f"./output/jatkasaari_network.csv"
    # Create output folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load dataset
    df = pd.read_csv(csv_path)

    plt.rcParams.update({
        "font.size": 13,
        "axes.labelsize": 16,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "legend.fontsize": 13
    })

    ################# Objective #################
    
    agg_df = df.groupby("budget")[["objective"]].first().reset_index()

    fig, ax = plt.subplots(figsize=(7, 3))
    ax.plot(agg_df["budget"], agg_df["objective"], marker="o", label="TTT")

    improv_perc = 100 - agg_df["objective"].iloc[-1] * 100 / agg_df["objective"].iloc[0]

    ax.set_xlabel("Budget [hours]")
    ax.set_ylabel("TTT [hours]")
    ax.set_xlim(left=0)
    ax.set_xlim(right=1e4)
    ax.grid(True, alpha=0.3)
    
    fig_path = os.path.join(output_dir, f"jatkasaari_network_TTT.png")
    fig.tight_layout()
    if saveFig:
        plt.savefig(fig_path, dpi=300)
        plt.close()
    else:
        plt.show()

    ################ qPT + qM #################

    OD_PT = {"3_13", "9_13"}
    OD_C = {"3_13", "9_13"}
    OD_V = {"23_24", "12_24", "14_16", "12_16"}

    agg_df_PT = (
        df[df["OD"].isin(OD_PT)]
        .groupby("budget")[["qpt"]]
        .sum()
        .reset_index()
    )

    agg_df_C= (
        df[df["OD"].isin(OD_C)]
        .groupby("budget")[["qc"]]
        .sum()
        .reset_index()
    )

    agg_df_M = (
        df[df["OD"].isin(OD_V)]
        .groupby("budget")[["qpt", "qc"]]
        .sum()
        .reset_index()
    )

    fig, ax = plt.subplots(figsize=(7, 3))
    ax.plot(agg_df_C["budget"], (agg_df_PT["qpt"]+agg_df_M["qpt"]) / (agg_df_PT["qpt"]+agg_df_C["qc"]+agg_df_M["qpt"]+agg_df_M["qc"]) * 100, marker="o", label="Total qpt (3_13 and 9_13)")
    ax.set_xlabel("Budget [hours]")
    ax.set_ylabel("Total usage of PT [%]")
    ax.set_xlim(left=0)
    ax.set_xlim(right=1e4)
    ax.grid(True, alpha=0.3)
    
    fig_path = os.path.join(output_dir, f"jatkasaari_network_PT_perc.png")
    fig.tight_layout()
    if saveFig:
        plt.savefig(fig_path, dpi=300)
        plt.close()
    else:
        plt.show()

    #print((agg_df_PT["qpt"]+agg_df_M["qpt"]) / (agg_df_PT["qpt"]+agg_df_C["qc"]+agg_df_M["qpt"]+agg_df_M["qc"]) * 100)

    ################### C/B ###################

    agg_df = (
        df.groupby("budget")[["objective"]]
        .first()
        .reset_index()
        .sort_values("budget")
    )
    #print(agg_df)
    baseline_obj = df[df["budget"] == 0]["objective"].iloc[0]
    agg_df["delta_obj"] = agg_df["objective"].shift(1) - agg_df["objective"]
    agg_df["delta_budget"] = agg_df["budget"] - agg_df["budget"].shift(1)

    agg_df["marginal_cb"] = agg_df["delta_obj"] / agg_df["delta_budget"]

    agg_df = (
        df.groupby("budget")[["objective"]]
        .first()
        .reset_index()
        .sort_values("budget")
    )

    agg_df = agg_df.sort_values("budget").reset_index(drop=True)
    agg_df["delta_obj"] = agg_df["objective"].shift(1) - agg_df["objective"]
    agg_df["delta_budget"] = agg_df["budget"] - agg_df["budget"].shift(1)

    agg_df["marginal_cb"] = agg_df["delta_obj"] / agg_df["delta_budget"]
    agg_df["pct_reduction"] = ((baseline_obj - agg_df["objective"]) / baseline_obj) * 100

    plot_df = agg_df.copy()

    # INTERPOLATE for continuous filling across full x-range
    x_full = np.linspace(plot_df["budget"].min(), plot_df["budget"].max(), 500)
    y_cb_interp = np.interp(x_full, plot_df["budget"], plot_df["marginal_cb"])

    # Plot with continuous interpolated filling
    fig, ax1 = plt.subplots(figsize=(7, 3))

    # Left: Total C/B ratio with CONTINUOUS filling
    ax1.plot(plot_df["budget"], plot_df["marginal_cb"], marker="o", linewidth=3, markersize=10, label="Total Benefit/Cost")
    ax1.axhline(y=1, color="red", linestyle="--", alpha=0.8, label="Break-even (Benefit = Cost)")

    # CONTINUOUS green fill: where interpolated C/B >= 1
    mask_green = y_cb_interp >= 1
    y_cb_interp[np.isposinf(y_cb_interp)] = 10
    ax1.fill_between(x_full, 1, y_cb_interp, where=mask_green, color="green", alpha=0.3, 
                    interpolate=True, label="Benefit > Cost", edgecolor="none")
    ax1.fill_between([0, 1000], 1, 1.5, color="green", alpha=0.3, edgecolor="none")

    # CONTINUOUS red fill: where interpolated C/B < 1  
    mask_red = y_cb_interp < 1
    ax1.fill_between(x_full, 1, y_cb_interp, where=mask_red, color="red", alpha=0.3, 
                    interpolate=True, label="Cost > Benefit")

    ax1.set_xlabel("Budget [hours]")
    ax1.set_ylabel("$\Delta$TTT / $\Delta$Budget")
    ax1.legend(loc="upper right", facecolor="white", framealpha=1)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(left=0)
    ax1.set_xlim(right=1e4)
    ax1.set_ylim(top=1.5)
    ax1.set_ylim(bottom=0)
    
    plt.tight_layout()
    fig_path = os.path.join(output_dir, f"jatkasaari_network_total_cost_benefit.png")
    fig.tight_layout()
    if saveFig:
        plt.savefig(fig_path, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()
    
    # Right: % reduction (unchanged)
    fig, ax2 = plt.subplots(figsize=(7, 3))
    ax2.plot(plot_df["budget"], plot_df["pct_reduction"], marker="s", color="blue")
    ax2.set_xlabel("Budget [hours]")
    ax2.set_ylabel("TTT reduction [%]")
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(left=0)
    ax2.set_xlim(right=1e4)
    ax2.set_ylim(bottom=0)
    
    plt.tight_layout()
    fig_path = os.path.join(output_dir, f"jatkasaari_network_TTT_reduction.png")
    fig.tight_layout()
    if saveFig:
        plt.savefig(fig_path, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()

    if saveFig:
        print("All plots generated successfully!")

if __name__ == "__main__":
    main()
