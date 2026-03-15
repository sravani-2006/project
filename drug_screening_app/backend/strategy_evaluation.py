from typing import Dict

import numpy as np
import pandas as pd

from fingerprints import generate_morgan_fingerprints, tanimoto_similarity_and_distance


def _pairwise_similarity_metrics(df: pd.DataFrame) -> tuple[float, int]:
    """
    Return average pairwise similarity and redundancy count (>0.8).
    """
    if df.empty or len(df) < 2:
        return 0.0, 0

    filtered_df, fp_matrix = generate_morgan_fingerprints(df)
    if len(filtered_df) < 2:
        return 0.0, 0

    similarity_matrix, _ = tanimoto_similarity_and_distance(fp_matrix)
    tri_upper = np.triu_indices_from(similarity_matrix, k=1)
    pair_values = similarity_matrix[tri_upper]

    if pair_values.size == 0:
        return 0.0, 0

    avg_similarity = float(pair_values.mean())
    redundancy_count = int((pair_values > 0.8).sum())
    return avg_similarity, redundancy_count


def evaluate_strategies(strategy_to_hits: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Evaluate strategies by average docking score, similarity, and redundancy.
    """
    rows = []

    for strategy_name, hits_df in strategy_to_hits.items():
        if hits_df is None or hits_df.empty:
            avg_docking = np.nan
            avg_similarity = 0.0
            redundancy_count = 0
        else:
            avg_docking = float(hits_df["docking_score"].astype(float).mean())
            avg_similarity, redundancy_count = _pairwise_similarity_metrics(hits_df)

        rows.append(
            {
                "strategy": strategy_name,
                "average_docking_score": avg_docking,
                "average_similarity": avg_similarity,
                "redundancy_count": redundancy_count,
            }
        )

    return pd.DataFrame(rows)
