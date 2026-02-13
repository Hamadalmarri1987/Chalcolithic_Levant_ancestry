import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def find_column(df: pd.DataFrame, candidates: list[str]) -> str:
    """Return the first matching column from candidates, case-insensitive."""
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    raise KeyError(f"None of these columns were found: {candidates}")


def coerce_numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def compute_dominance_ranges(
    df: pd.DataFrame,
    p_col: str,
    prop_cols: dict[str, str],
    se_cols: dict[str, str],
    p_min: float = 0.60,
    p_max: float = 0.999,
    se_max: float = 0.08,
) -> tuple[pd.DataFrame, pd.DataFrame]:

    numeric_cols = [p_col] + list(prop_cols.values()) + list(se_cols.values())
    numeric_cols = [c for c in numeric_cols if c in df.columns]
    df = df.copy()
    df = coerce_numeric(df, numeric_cols)

    # ---------------------------------
    # NO p-VALUE FILTERING (ONLY CHANGE)
    # ---------------------------------
    f = df.copy()

    se_present = f[list(se_cols.values())].copy()
    any_bad_se = (se_present > se_max).any(axis=1)
    f = f[~any_bad_se].copy()

    out_rows = []
    for comp, prop_c in prop_cols.items():
        se_c = se_cols.get(comp, None)

        sub = f[~f[prop_c].isna()].copy()
        if se_c is not None and se_c in sub.columns:
            sub = sub[~sub[se_c].isna()]

        if sub.empty:
            out_rows.append([comp, np.nan, np.nan, 0])
            continue

        out_rows.append([
            comp,
            float(sub[prop_c].min()),
            float(sub[prop_c].max()),
            int(len(sub))
        ])

    ranges_df = pd.DataFrame(
        out_rows,
        columns=["component", "min_prop", "max_prop", "n_models_used"]
    )
    return ranges_df, f


COLOR_MAP = {
    "Levant": "#2ca02c",
    "Anatolia": "#1f77b4",
    "Mesopotamia": "#9467bd",
    "Iran/Zagros": "#d62728",
    "ZagrosMeso": "#ff7f0e",
}


def plot_dominance_ranges(
    ranges_df: pd.DataFrame,
    order: list[str],
    title: str,
    out_png: str,
):
    ranges_df = ranges_df.set_index("component").reindex(order).reset_index()
    plot_df = ranges_df.dropna(subset=["min_prop", "max_prop"]).copy()

    fig, ax = plt.subplots(figsize=(10, 6))

    y = np.arange(len(plot_df))
    left = plot_df["min_prop"].values
    width = (plot_df["max_prop"] - plot_df["min_prop"]).values
    colors = [COLOR_MAP[c] for c in plot_df["component"]]

    ax.barh(y, width, left=left, color=colors)

    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["component"], fontsize=9)

    ax.set_xlabel("Ancestry proportion", fontsize=14, fontweight="bold")
    ax.set_title(title, fontsize=16, fontweight="bold", pad=5)

    ax.tick_params(axis="x", labelsize=13)
    ax.tick_params(axis="y", labelsize=14)

    ax.invert_yaxis()

    for i, row in plot_df.reset_index(drop=True).iterrows():
        a = int(row["min_prop"] * 100)
        b = int(row["max_prop"] * 100)
        n = int(row["n_models_used"])

        ax.text(
            row["max_prop"] + 0.01,
            i,
            f"{a}–{b}% (n={n})",
            va="center",
            fontsize=11
        )

    ax.set_xlim(0, 0.8)
    ax.grid(True, axis="x", alpha=0.25)

    plt.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main():
    excel_path = r"C:\\Users\\WAQAS\\Downloads\\graphs_for_research\\Copy of S1_qpAdm.xlsx"
    sheet_name = "3way_passing_models"

    p_min = 0.60
    p_max = 0.999
    se_max = 0.08

    out_png = "dominance_ranges_allp_se008.png"
    out_csv = "filtered_models_allp_se008.csv"

    if not os.path.exists(excel_path):
        raise FileNotFoundError(f"Excel file not found: {excel_path}")

    df = pd.read_excel(excel_path, sheet_name=sheet_name)

    p_col = find_column(df, ["p", "p_value", "pval", "p-value"])

    prop_cols = {
        "Levant": find_column(df, ["Levant_prop", "Levant"]),
        "Anatolia": find_column(df, ["Anatolia_prop", "Anatolia"]),
        "Mesopotamia": find_column(df, ["Mesopotamia_prop", "Mesopotamia"]),
        "Iran/Zagros": find_column(df, ["Iran_prop", "Iran"]),
        "ZagrosMeso": find_column(df, ["ZagrosMeso_prop", "ZagrosMeso"]),
    }

    se_cols = {
        "Levant": find_column(df, ["Levant_se"]),
        "Anatolia": find_column(df, ["Anatolia_se"]),
        "Mesopotamia": find_column(df, ["Mesopotamia_se"]),
        "Iran/Zagros": find_column(df, ["Iran_se"]),
        "ZagrosMeso": find_column(df, ["ZagrosMeso_se"]),
    }

    ranges_df, filtered_df = compute_dominance_ranges(
        df=df,
        p_col=p_col,
        prop_cols=prop_cols,
        se_cols=se_cols,
        p_min=p_min,
        p_max=p_max,
        se_max=se_max,
    )

    filtered_df.to_csv(out_csv, index=False)

    order = ["Levant", "Anatolia", "Mesopotamia", "Iran/Zagros", "ZagrosMeso"]
    title = f"Dominant ancestry ranges (all SE ≤ {se_max})"


    plot_dominance_ranges(
        ranges_df=ranges_df,
        order=order,
        title=title,
        out_png=out_png,
    )

    print("Dominance ranges (from filtered models):")
    for _, r in ranges_df.dropna(subset=["min_prop", "max_prop"]).iterrows():
        print(
            f"- {r['component']}: {r['min_prop']*100:.1f}–{r['max_prop']*100:.1f}% "
            f"(n={int(r['n_models_used'])})"
        )

    print(f"\nSaved figure: {out_png}")
    print(f"Saved filtered table: {out_csv}")


if __name__ == "__main__":
    main()
