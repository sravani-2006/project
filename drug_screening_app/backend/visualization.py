import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


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
