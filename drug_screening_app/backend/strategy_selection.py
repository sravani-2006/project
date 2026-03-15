import numpy as np
import pandas as pd

from fingerprints import generate_morgan_fingerprints, tanimoto_similarity_and_distance


def select_score_only_hits(docking_df: pd.DataFrame, k: int = 20) -> pd.DataFrame:
    """Baseline selection using only docking score (lower is better)."""
    if docking_df.empty:
        return pd.DataFrame(columns=["ligand_id", "smiles", "docking_score"])

    out = (
        docking_df.sort_values("docking_score", ascending=True)
        .head(k)
        .reset_index(drop=True)
    )
    return out[["ligand_id", "smiles", "docking_score"]]


def select_greedy_diversity_hits(
    docking_df: pd.DataFrame,
    k: int = 20,
    radius: int = 2,
    n_bits: int = 1024,
) -> pd.DataFrame:
    """
    Greedy max-min diversity selection using Tanimoto distance.

    Steps:
    - First pick the best docking score molecule
    - Then iteratively add the molecule with maximum minimum distance
      from the currently selected set.
    """
    if docking_df.empty:
        return pd.DataFrame(columns=["ligand_id", "smiles", "docking_score"])

    filtered_df, fp_matrix = generate_morgan_fingerprints(docking_df, radius=radius, n_bits=n_bits)
    if filtered_df.empty:
        return pd.DataFrame(columns=["ligand_id", "smiles", "docking_score"])

    _, distance_matrix = tanimoto_similarity_and_distance(fp_matrix)
    n = len(filtered_df)
    target_k = min(k, n)

    first_idx = int(filtered_df["docking_score"].astype(float).idxmin())
    selected = [first_idx]

    while len(selected) < target_k:
        remaining = [i for i in range(n) if i not in selected]

        best_idx = None
        best_min_dist = -1.0
        best_score = float("inf")

        for i in remaining:
            min_dist = float(distance_matrix[i, selected].min())
            score_i = float(filtered_df.iloc[i]["docking_score"])

            if min_dist > best_min_dist or (np.isclose(min_dist, best_min_dist) and score_i < best_score):
                best_min_dist = min_dist
                best_score = score_i
                best_idx = i

        if best_idx is None:
            break

        selected.append(best_idx)

    out = filtered_df.iloc[selected].copy()
    out = out.sort_values("docking_score", ascending=True).reset_index(drop=True)
    return out[["ligand_id", "smiles", "docking_score"]]


def select_multiobjective_hits(
    docking_df: pd.DataFrame,
    k: int = 20,
    radius: int = 2,
    n_bits: int = 1024,
) -> pd.DataFrame:
    """
    Multi-objective selection balancing docking and diversity.

    final_score = 0.7 * normalized_binding + 0.3 * diversity
    where normalized_binding maps lower docking score to higher value.
    """
    if docking_df.empty:
        return pd.DataFrame(columns=["ligand_id", "smiles", "docking_score", "diversity_score", "final_score"])

    filtered_df, fp_matrix = generate_morgan_fingerprints(docking_df, radius=radius, n_bits=n_bits)
    if filtered_df.empty:
        return pd.DataFrame(columns=["ligand_id", "smiles", "docking_score", "diversity_score", "final_score"])

    _, distance_matrix = tanimoto_similarity_and_distance(fp_matrix)

    scores = filtered_df["docking_score"].astype(float).to_numpy()
    s_min = float(scores.min())
    s_max = float(scores.max())

    if np.isclose(s_max, s_min):
        normalized_binding = np.ones_like(scores, dtype=float)
    else:
        normalized_binding = (s_max - scores) / (s_max - s_min)

    n = len(filtered_df)
    if n <= 1:
        diversity_scores = np.zeros(n, dtype=float)
    else:
        diversity_scores = distance_matrix.sum(axis=1) / (n - 1)

    final_scores = 0.7 * normalized_binding + 0.3 * diversity_scores

    out = filtered_df.copy()
    out["diversity_score"] = diversity_scores
    out["final_score"] = final_scores

    out = out.sort_values(["final_score", "docking_score"], ascending=[False, True]).head(k).reset_index(drop=True)

    return out[["ligand_id", "smiles", "docking_score", "diversity_score", "final_score"]]
