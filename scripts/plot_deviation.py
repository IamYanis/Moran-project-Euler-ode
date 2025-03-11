import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data
file_path = "data/moran_simulation_results.csv"
df = pd.read_csv(file_path)

# Filter deviation data
df_dev = df[df["Type"] == "Deviation"].copy()

# Number of simulations (as used in main.c)
N_SIM = 1000  
Z_95 = 1.96  # 95% confidence interval Z-score

# Compute 95% CI for error bars
df_dev["CI_95"] = (df_dev["std_sup_dev"] / np.sqrt(N_SIM)) * Z_95

# Plot mean supremum deviation vs N with proper error bars
plt.figure(figsize=(8, 5))
for r in sorted(df_dev["r"].unique()):
    subset = df_dev[df_dev["r"] == r]
    plt.errorbar(subset["N"], subset["Value"], yerr=subset["CI_95"], label=f'r={r:.2f}', capsize=2)

plt.xlabel("Population size N")
plt.ylabel("Mean Supremum Deviation")
plt.title("Mean Supremum Deviation vs N with 95% Confidence Intervals")
plt.legend()
plt.grid()
plt.show()
