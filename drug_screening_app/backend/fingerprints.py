import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.metrics import pairwise_distances


def generate_morgan_fingerprints(hit_df: pd.DataFrame, radius: int = 2, n_bits: int = 1024):
    """
    Generate Morgan fingerprints for molecules.

    Returns:
    - filtered DataFrame (invalid SMILES removed)
    - fingerprint matrix with shape [n_molecules, n_bits]
    """
    vectors = []
    valid_rows = []

    print("[Fingerprints] Generating Morgan fingerprints...")

    for _, row in hit_df.iterrows():
        smiles = str(row["smiles"])
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        arr = np.zeros((n_bits,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)

        vectors.append(arr)
        valid_rows.append(row.to_dict())

    filtered_df = pd.DataFrame(valid_rows)
    matrix = np.array(vectors, dtype=np.uint8)

    print(f"[Fingerprints] Matrix shape: {matrix.shape}")
    return filtered_df, matrix


def tanimoto_similarity_and_distance(fingerprint_matrix: np.ndarray):
    """
    Compute Tanimoto similarity and distance matrices from binary fingerprints.

    For binary vectors, Jaccard distance is equivalent to Tanimoto distance.
    """
    if len(fingerprint_matrix) == 0:
        return np.empty((0, 0)), np.empty((0, 0))

    # sklearn expects boolean-like data for Jaccard.
    bool_matrix = fingerprint_matrix.astype(bool)
    distance_matrix = pairwise_distances(bool_matrix, metric="jaccard")
    similarity_matrix = 1.0 - distance_matrix

    print(f"[Similarity] Distance matrix shape: {distance_matrix.shape}")
    return similarity_matrix, distance_matrix
