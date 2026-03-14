import pandas as pd
from sklearn.cluster import DBSCAN


def cluster_hits(
    hit_df: pd.DataFrame,
    fingerprint_matrix,
    eps: float = 0.7,
    min_samples: int = 2,
) -> pd.DataFrame:
    """
    Cluster molecules using DBSCAN with Jaccard distance on binary fingerprints.
    """
    print("[Clustering] Running DBSCAN clustering...")
    print(f"[Clustering] Number of molecules loaded: {len(hit_df)}")
    print(f"[Clustering] Number of fingerprints generated: {len(fingerprint_matrix)}")
    print(f"[Clustering] Parameters -> eps={eps}, min_samples={min_samples}, metric=jaccard")

    if len(hit_df) < 2:
        print("[Clustering] WARNING: Not enough molecules to cluster (need at least 2).")
        out = hit_df.copy()
        out["cluster"] = -1
        return out

    model = DBSCAN(eps=eps, min_samples=min_samples, metric="jaccard")
    labels = model.fit_predict(fingerprint_matrix)

    clustered = hit_df.copy()
    clustered["cluster"] = labels

    distribution = pd.Series(labels).value_counts().sort_index().to_dict()
    print(f"[Clustering] Cluster labels distribution: {distribution}")

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = int((labels == -1).sum())

    print(f"[Clustering] Clusters found: {n_clusters}")
    print(f"[Clustering] Noise molecules: {n_noise}")
    if n_clusters == 0:
        print(
            "[Clustering] WARNING: No clusters found. Try increasing eps (e.g. 0.8-0.9) "
            "or reducing min_samples."
        )
    return clustered


def cluster_hits_from_distance(
    hit_df: pd.DataFrame,
    distance_matrix,
    eps: float = 0.7,
    min_samples: int = 2,
) -> pd.DataFrame:
    """
    Cluster using a precomputed Tanimoto/Jaccard distance matrix.
    """
    print("[Clustering] Running DBSCAN on precomputed distance matrix...")
    print(f"[Clustering] Number of molecules loaded: {len(hit_df)}")

    if len(hit_df) < 2:
        out = hit_df.copy()
        out["cluster"] = -1
        print("[Clustering] WARNING: Not enough molecules to cluster.")
        return out

    model = DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed")
    labels = model.fit_predict(distance_matrix)

    clustered = hit_df.copy()
    clustered["cluster"] = labels

    distribution = pd.Series(labels).value_counts().sort_index().to_dict()
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = int((labels == -1).sum())

    print(f"[Clustering] Cluster labels distribution: {distribution}")
    print(f"[Clustering] Number of clusters found: {n_clusters}")
    print(f"[Clustering] Noise molecules: {n_noise}")

    if n_clusters == 0:
        print(
            "[Clustering] WARNING: No clusters found. Suggestion: increase eps or reduce min_samples."
        )

    return clustered
