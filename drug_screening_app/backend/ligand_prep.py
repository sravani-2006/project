import os
import subprocess
import tempfile
from typing import Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


REQUIRED_COLUMNS = {"ligand_id", "smiles"}


def _smiles_to_3d_mol(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    embed_status = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if embed_status != 0:
        return None

    AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    return mol


def _convert_mol_to_pdbqt(mol, output_path: str) -> bool:
    # RDKit handles chemistry; Open Babel converts MOL to PDBQT format for Vina.
    with tempfile.NamedTemporaryFile(suffix=".mol", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        Chem.MolToMolFile(mol, tmp_path)
        command = ["obabel", tmp_path, "-O", output_path]
        subprocess.run(command, check=True, capture_output=True, text=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"[Ligand Prep] Failed PDBQT conversion for {output_path}: {exc}")
        return False
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def prepare_ligands(ligands_csv_path: str, output_dir: str, max_ligands: Optional[int] = None) -> pd.DataFrame:
    """
    Convert ligand SMILES into 3D structures and PDBQT files.

    Returns a DataFrame with ligand_id, smiles, and pdbqt_path for valid ligands.
    """
    os.makedirs(output_dir, exist_ok=True)

    print("[Ligand Prep] Loading ligand CSV...")
    ligands_df = pd.read_csv(ligands_csv_path)

    missing = REQUIRED_COLUMNS - set(ligands_df.columns)
    if missing:
        raise ValueError(f"Ligand CSV is missing columns: {sorted(missing)}")

    if max_ligands is not None:
        ligands_df = ligands_df.head(max_ligands)

    prepared_rows = []
    invalid_smiles = 0
    conversion_failures = 0

    total = len(ligands_df)
    print(f"[Ligand Prep] Preparing {total} ligands...")

    for idx, row in ligands_df.iterrows():
        ligand_id = str(row["ligand_id"])
        smiles = str(row["smiles"])

        mol = _smiles_to_3d_mol(smiles)
        if mol is None:
            invalid_smiles += 1
            continue

        pdbqt_path = os.path.join(output_dir, f"{ligand_id}.pdbqt")
        ok = _convert_mol_to_pdbqt(mol, pdbqt_path)
        if not ok:
            conversion_failures += 1
            continue

        prepared_rows.append({
            "ligand_id": ligand_id,
            "smiles": smiles,
            "pdbqt_path": pdbqt_path,
        })

        if (idx + 1) % 500 == 0:
            print(f"[Ligand Prep] Processed {idx + 1}/{total} ligands...")

    prepared_df = pd.DataFrame(prepared_rows)

    print(
        "[Ligand Prep] Completed. "
        f"Prepared: {len(prepared_df)}, Invalid SMILES: {invalid_smiles}, Conversion failures: {conversion_failures}"
    )

    prepared_df.to_csv(os.path.join(output_dir, "prepared_ligands.csv"), index=False)
    return prepared_df
