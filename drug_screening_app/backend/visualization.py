import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from fingerprints import generate_morgan_fingerprints


def save_docking_score_distribution(docking_df: pd.DataFrame, output_path: str) -> str:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    plt.figure(figsize=(8, 4))
    plt.hist(docking_df["docking_score"], bins=20, color="#d44d2c", edgecolor="white")
    plt.xlabel("Docking Score")
    plt.ylabel("Count")
    plt.title("Docking Score Distribution")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    return output_path


def save_cluster_distribution(clustered_df: pd.DataFrame, output_path: str) -> str:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    counts = clustered_df[clustered_df["cluster"] != -1]["cluster"].value_counts().sort_index()

    plt.figure(figsize=(10, 4))
    if len(counts) > 0:
        plt.bar(counts.index.astype(str), counts.values, color="#f09b3f")
    plt.xlabel("Cluster")
    plt.ylabel("Molecule Count")
    plt.title("Cluster Size Distribution")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    return output_path


def save_chemical_space_pca(fingerprint_matrix: np.ndarray, labels, output_path: str) -> str:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Handle tiny inputs gracefully.
    if len(fingerprint_matrix) < 2:
        plt.figure(figsize=(7, 5))
        plt.text(0.5, 0.5, "Not enough molecules for PCA", ha="center", va="center")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()
        return output_path

    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(fingerprint_matrix)

    plt.figure(figsize=(7, 5))
    scatter = plt.scatter(coords[:, 0], coords[:, 1], c=labels, cmap="tab20", s=35, alpha=0.85)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Chemical Space (PCA)")
    plt.colorbar(scatter, label="Cluster")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    return output_path


def save_strategy_metric_comparison(
    evaluation_df: pd.DataFrame,
    metric_column: str,
    title: str,
    y_label: str,
    output_path: str,
) -> str:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    plot_df = evaluation_df.copy()
    if plot_df.empty or metric_column not in plot_df.columns:
        plt.figure(figsize=(8, 4))
        plt.text(0.5, 0.5, "No data for comparison", ha="center", va="center")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()
        return output_path

    labels = plot_df["strategy"].astype(str).tolist()
    values = plot_df[metric_column].astype(float).tolist()

    plt.figure(figsize=(9, 4.8))
    bars = plt.bar(labels, values, color=["#d44d2c", "#f09b3f", "#5a9bd5", "#2f7d58"][: len(labels)])
    for bar, value in zip(bars, values):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{value:.3f}", ha="center", va="bottom", fontsize=9)

    plt.title(title)
    plt.xlabel("Strategy")
    plt.ylabel(y_label)
    plt.xticks(rotation=15)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    return output_path


def save_strategy_chemical_space_pca(strategy_to_hits: dict[str, pd.DataFrame], output_path: str) -> str:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    rows = []
    for strategy, df in strategy_to_hits.items():
        if df is None or df.empty:
            continue
        mini = df[["ligand_id", "smiles"]].copy()
        mini["strategy"] = strategy
        rows.append(mini)

    if not rows:
        plt.figure(figsize=(7, 5))
        plt.text(0.5, 0.5, "No selected molecules for PCA", ha="center", va="center")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()
        return output_path

    combined = pd.concat(rows, ignore_index=True).drop_duplicates(subset=["ligand_id", "smiles"])
    fp_df, fp_matrix = generate_morgan_fingerprints(combined)

    if len(fp_df) < 2:
        plt.figure(figsize=(7, 5))
        plt.text(0.5, 0.5, "Not enough molecules for PCA", ha="center", va="center")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()
        return output_path

    plot_df = fp_df.copy()
    if "strategy" not in plot_df.columns:
        plot_df = plot_df.merge(combined[["ligand_id", "smiles", "strategy"]], on=["ligand_id", "smiles"], how="left")

    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(fp_matrix)

    plt.figure(figsize=(7.2, 5.4))
    for strategy in plot_df["strategy"].dropna().unique():
        idx = plot_df["strategy"] == strategy
        plt.scatter(coords[idx, 0], coords[idx, 1], s=35, alpha=0.85, label=strategy)

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Chemical Space PCA (Selected Hits by Strategy)")
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    return output_path
