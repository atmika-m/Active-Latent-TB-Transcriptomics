############################################################
# Random Forest-Based Biomarker Prioritization
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
# This script compares multiple regression-based machine
# learning models using differential gene expression
# statistics and prioritizes candidate biomarkers using
# Random Forest prediction scores.
############################################################

# ==========================================================
# Import Required Libraries
# ==========================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import mean_squared_error, r2_score


# ==========================================================
# Load Dataset
# ==========================================================

# Load merged annotated DEG table
data = pd.read_csv("../results/annotated/AMC.csv")


# ==========================================================
# Prepare Features and Target Variable
# ==========================================================

# Target variable
target = "logFC"

# Predictor variables
X = data[
    [
        "AveExpr",
        "t",
        "P.Value",
        "B",
        "adj.P.Val"
    ]
]

y = data[target]


# ==========================================================
# Split Dataset
# ==========================================================

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.20,
    random_state=42
)


# ==========================================================
# Train Machine Learning Models
# ==========================================================

# Random Forest
rf = RandomForestRegressor(
    n_estimators=100,
    random_state=42
)

rf.fit(X_train, y_train)

y_pred_rf = rf.predict(X_test)


# Gradient Boosting
gb = GradientBoostingRegressor(
    n_estimators=100,
    learning_rate=0.1,
    max_depth=3,
    random_state=42
)

gb.fit(X_train, y_train)

y_pred_gb = gb.predict(X_test)


# Decision Tree
dt = DecisionTreeRegressor(
    max_depth=5,
    random_state=42
)

dt.fit(X_train, y_train)

y_pred_dt = dt.predict(X_test)


# ==========================================================
# Evaluate Model Performance
# ==========================================================

print("\nDecision Tree Results")
print("----------------------")
print("MSE:", mean_squared_error(y_test, y_pred_dt))
print("R²:", r2_score(y_test, y_pred_dt))

print("\nRandom Forest Results")
print("----------------------")
print("MSE:", mean_squared_error(y_test, y_pred_rf))
print("R²:", r2_score(y_test, y_pred_rf))

print("\nGradient Boosting Results")
print("--------------------------")
print("MSE:", mean_squared_error(y_test, y_pred_gb))
print("R²:", r2_score(y_test, y_pred_gb))


# ==========================================================
# Compare Model Performance
# ==========================================================

results = pd.DataFrame({

    "Model": [

        "Decision Tree",

        "Random Forest",

        "Gradient Boosting"

    ],

    "R2 Score": [

        r2_score(y_test, y_pred_dt),

        r2_score(y_test, y_pred_rf),

        r2_score(y_test, y_pred_gb)

    ]

})

print("\nModel Performance Summary")
print(results)


# ==========================================================
# Visualize Prediction Performance
# ==========================================================

plt.figure(figsize=(7,5))

plt.scatter(
    y_test,
    y_pred_rf,
    label="Random Forest",
    color="red"
)

plt.scatter(
    y_test,
    y_pred_gb,
    label="Gradient Boosting",
    color="blue"
)

plt.scatter(
    y_test,
    y_pred_dt,
    label="Decision Tree",
    color="deeppink"
)

plt.plot(

    [min(y_test), max(y_test)],

    [min(y_test), max(y_test)],

    "k--"

)

plt.xlabel("Observed logFC")

plt.ylabel("Predicted logFC")

plt.title("Regression Model Predictions")

plt.legend()

plt.tight_layout()

plt.show()


# ==========================================================
# Compare R² Scores
# ==========================================================

models = {

    "Random Forest": rf,

    "Decision Tree": dt,

    "Gradient Boosting": gb

}

results = {}

for name, model in models.items():

    prediction = model.predict(X_test)

    results[name] = {

        "MSE": mean_squared_error(y_test, prediction),

        "R2": r2_score(y_test, prediction)

    }

    print(
        f"{name}: "
        f"MSE = {results[name]['MSE']:.4f}, "
        f"R² = {results[name]['R2']:.4f}"
    )


plt.figure(figsize=(8,5))

plt.barh(

    list(results.keys()),

    [results[m]["R2"] for m in results],

    color="steelblue"

)

plt.xlabel("R² Score")

plt.title("Machine Learning Model Comparison")

plt.tight_layout()

plt.show()


# ==========================================================
# Random Forest Biomarker Prioritization
# ==========================================================

# Predict Random Forest scores
data["RF_score"] = rf.predict(

    data[
        [
            "AveExpr",
            "t",
            "P.Value",
            "B",
            "adj.P.Val"
        ]
    ]

)

# Rank genes
top50 = (

    data

    .sort_values(

        by="RF_score",

        ascending=False

    )

    .head(50)

)

print("\nTop 50 Candidate Biomarkers")

print(top50[["GeneSymbol", "RF_score"]])


# ==========================================================
# Save Results
# ==========================================================

top50[["GeneSymbol", "RF_score"]].to_csv(

    "../results/machine_learning/RF_Top50_Genes.csv",

    index=False

)

print("\nRandom Forest biomarker prioritization completed successfully.")
