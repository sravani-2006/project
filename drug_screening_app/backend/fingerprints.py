import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


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
