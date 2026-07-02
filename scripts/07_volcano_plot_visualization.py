############################################################
# Volcano Plot Visualization
#
# Project:
# Comparative Transcriptomic Profiling of Active and
# Latent Tuberculosis to Identify Differential
# Host-Gene Expression Signatures
#
# Author:
# Mrunalline Atmika
#
# Description:
# This script generates volcano plots for each
# transcriptomic dataset using differential gene
# expression results obtained from limma analysis.
#
# Datasets:
# - GSE19491
# - GSE28623
# - GSE25534
############################################################

# ==========================================================
# Import Required Libraries
# ==========================================================

import numpy as np
import matplotlib.pyplot as plt


# ==========================================================
# Store Differential Expression Results
# ==========================================================

datasets = {
    "GSE19491": deg1,
    "GSE28623": deg3,
    "GSE25534": deg4
}


# ==========================================================
# Generate Volcano Plot for Each Dataset
# ==========================================================

for dataset_name, deg_data in datasets.items():

    # Create a copy to avoid modifying the original dataframe
    deg = deg_data.copy()

    # ------------------------------------------------------
    # Classify Differentially Expressed Genes
    # ------------------------------------------------------

    deg["Significance"] = "Not Significant"

    deg.loc[
        (deg["adj.P.Val"] < 0.05) &
        (deg["logFC"] > 1),
        "Significance"
    ] = "Upregulated"

    deg.loc[
        (deg["adj.P.Val"] < 0.05) &
        (deg["logFC"] < -1),
        "Significance"
    ] = "Downregulated"

    # ------------------------------------------------------
    # Convert adjusted p-values for plotting
    # ------------------------------------------------------

    deg["adj.P.Val"] = deg["adj.P.Val"].replace(
        0,
        np.nextafter(0, 1)
    )

    deg["neg_log10_padj"] = -np.log10(deg["adj.P.Val"])

    # ------------------------------------------------------
    # Create Volcano Plot
    # ------------------------------------------------------

    plt.figure(figsize=(8, 6))

    colors = {
        "Upregulated": "red",
        "Downregulated": "blue",
        "Not Significant": "lightgrey"
    }

    for category, color in colors.items():

        subset = deg[deg["Significance"] == category]

        plt.scatter(
            subset["logFC"],
            subset["neg_log10_padj"],
            c=color,
            label=category,
            alpha=0.6,
            s=15
        )

    # ------------------------------------------------------
    # Add Threshold Lines
    # ------------------------------------------------------

    plt.axvline(
        x=1,
        color="black",
        linestyle="--"
    )

    plt.axvline(
        x=-1,
        color="black",
        linestyle="--"
    )

    plt.axhline(
        y=-np.log10(0.05),
        color="black",
        linestyle="--"
    )

    # ------------------------------------------------------
    # Add Labels and Title
    # ------------------------------------------------------

    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10 Adjusted P-value")
    plt.title(f"Volcano Plot - {dataset_name}")

    plt.legend()

    plt.tight_layout()

    # ------------------------------------------------------
    # Save Figure
    # ------------------------------------------------------

    plt.savefig(
        f"../figures/Volcano_{dataset_name}.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()

print("Volcano plots generated successfully.")
