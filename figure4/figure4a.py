# -*- coding: utf-8 -*-
# qpAdm 3-way plot
# Clean publication-ready version

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# CONFIG
# -----------------------------
XLSX_PATH = r"path/to/qpadm_excel"
SHEET_NAME = "3way_passing_models"

KEEP_OPAQUE_1BASED = {2, 5, 6}

COLORS = {
    "L": {"dark": "#2ca02c", "light": "#ACFDAC"},
    "A": {"dark": "#1f77b4", "light": "#a0d6fc"},
    "M": {"dark": "#9467bd", "light": "#e7cbff"},
}

COMPONENTS = [
    ("L", ["Levant_prop", "Levant", "L_prop"], ["Levant_se", "L_se"]),
    ("A", ["Anatolia_prop", "Anatolia", "A_prop"], ["Anatolia_se", "A_se"]),
    ("M", ["Mesopotamia_prop", "Meso_prop", "E_prop"],
          ["Mesopotamia_se", "Meso_se", "E_se"]),
]

# -----------------------------
# Helpers
# -----------------------------
def pick_column(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"Missing columns from {candidates}")

def fmt_p(p):
    return f"{p:.4f}"

def fmt_se(x):
    return f"{int(x*1000)/1000:.3f}"  # truncate, no rounding

# -----------------------------
# Load & clean
# -----------------------------
df = pd.read_excel(XLSX_PATH, sheet_name=SHEET_NAME)

col_p = pick_column(df, ["p", "pval", "p_value", "p-value"])

resolved = {}
for key, prop_candidates, se_candidates in COMPONENTS:
    resolved[key] = {
        "prop": pick_column(df, prop_candidates),
        "se": pick_column(df, se_candidates),
    }

# Convert to numeric
for col in [col_p] + [v for comp in resolved.values() for v in comp.values()]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df = df.dropna(subset=[col_p])
df = df[df[col_p] >= 0.99]

df = df.dropna(subset=[v for comp in resolved.values() for v in comp.values()])

df = df.sort_values(col_p, ascending=False).reset_index(drop=True)
df["prominent"] = [(i + 1) in KEEP_OPAQUE_1BASED for i in range(len(df))]

# -----------------------------
# Plot
# -----------------------------
n = len(df)
x = np.arange(n)
width = 0.25
offsets = {"L": -width, "A": 0, "M": width}

fig, ax = plt.subplots(figsize=(14, 6))

for i, row in df.iterrows():
    dom = row["prominent"]
    bar_alpha = 1.0 if dom else 0.35
    err_alpha = 0.9 if dom else 0.35

    error_kw = dict(ecolor="gray", elinewidth=1.3,
                    capthick=1.3, alpha=err_alpha)

    for key in ["L", "A", "M"]:
        prop = row[resolved[key]["prop"]]
        se = row[resolved[key]["se"]]

        color = COLORS[key]["dark"] if dom else COLORS[key]["light"]

        xpos = x[i] + offsets[key]

        ax.bar(xpos, prop, width,
               yerr=se,
               color=color,
               alpha=bar_alpha,
               capsize=6,
               error_kw=error_kw)

        # Percent labels (automatic rounding)
        if dom:
            pct = int(round(prop * 100))
            ax.text(xpos,
                    prop + se + 0.01,
                    f"{pct}%",
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    fontweight="bold")

# -----------------------------
# X-axis labels
# -----------------------------
ax.set_xticks(x)
ax.set_xticklabels([
    f"p={fmt_p(row[col_p])}\n"
    f"L({fmt_se(row[resolved['L']['se']])})\n"
    f"A({fmt_se(row[resolved['A']['se']])})\n"
    f"M({fmt_se(row[resolved['M']['se']])})"
    for _, row in df.iterrows()
], fontsize=9)

# -----------------------------
# Cosmetics
# -----------------------------
ax.set_title(
    "All qpAdm models with p ≥ 0.99\n"
    "Optimal models bars: lower SE (≤ 0.08)\n"
    "Blurred bars: higher SE (> 0.08)",
    fontsize=10, fontweight="bold", pad=5
)

ax.set_ylabel("Ancestry proportion", fontsize=10, fontweight="bold")
ax.set_ylim(0, 0.6)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="y", labelsize=8)

# Legend inside (top-left)
ax.legend(
    handles=[
        plt.Rectangle((0, 0), 1, 1, color=COLORS["L"]["dark"]),
        plt.Rectangle((0, 0), 1, 1, color=COLORS["A"]["dark"]),
        plt.Rectangle((0, 0), 1, 1, color=COLORS["M"]["dark"]),
    ],
    labels=["Levant", "Anatolia", "Mesopotamia"],
    loc="upper left",
    frameon=False,
    fontsize=10
)

fig.tight_layout()
plt.show()

