import pandas as pd


def filter_top_by_docking_score(docking_df: pd.DataFrame, top_n: int = 8000) -> pd.DataFrame:
    """Select top molecules with the lowest docking scores."""
    print(f"[Filtering] Selecting top {top_n} molecules by docking score...")
    top_hits = docking_df.sort_values("docking_score", ascending=True).head(top_n).reset_index(drop=True)
    return top_hits


def select_cluster_representatives(clustered_df: pd.DataFrame) -> pd.DataFrame:
    """
    Select one representative per cluster (excluding noise cluster -1).

    Representative rule:
    - Molecule with the lowest docking score within each cluster.
    """
    print("[Hit Selection] Selecting representative molecules from each cluster...")

    non_noise = clustered_df[clustered_df["cluster"] != -1].copy()
    if non_noise.empty:
        return pd.DataFrame(columns=["ligand_id", "smiles", "docking_score", "cluster"])

    reps = (
        non_noise.sort_values("docking_score", ascending=True)
        .groupby("cluster", as_index=False)
        .first()[["ligand_id", "smiles", "docking_score", "cluster"]]
        .sort_values("docking_score", ascending=True)
        .reset_index(drop=True)
    )

    print(f"[Hit Selection] Final representative hits: {len(reps)}")
    return reps
