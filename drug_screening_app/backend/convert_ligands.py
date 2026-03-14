"""
Convert ligand files to the required CSV format (ligand_id, smiles).

Supported input formats:
  - PDBQT file (.pdbqt)           -> extracts SMILES from REMARK lines if present,
                                     or converts structures via RDKit
  - Multi-molecule SDF (.sdf)     -> reads directly with RDKit
  - CSV with wrong column names   -> lets you remap columns interactively
  - SMILES text file (.smi/.txt)  -> one SMILES per line, optional name column

Usage:
  python convert_ligands.py --input <your_file> --output ligands.csv
"""

import argparse
import os
import sys

import pandas as pd


def _from_pdbqt(path: str) -> pd.DataFrame:
    """
    Extract molecules from a multi-molecule PDBQT file.
    Looks for REMARK SMILES lines written by Open Babel.
    Falls back to reading molecule blocks and converting via RDKit.
    """
    from rdkit import Chem

    rows = []
    current_smiles = None
    current_name = None
    mol_block_lines = []
    mol_index = 0

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line_stripped = line.strip()

            # Open Babel writes: REMARK  Name = <name>
            if line_stripped.startswith("REMARK  Name ="):
                current_name = line_stripped.split("=", 1)[1].strip()

            # Open Babel writes: REMARK SMILES <smiles> <name>
            if line_stripped.startswith("REMARK SMILES "):
                parts = line_stripped.split(None, 3)
                if len(parts) >= 3:
                    current_smiles = parts[2]

            mol_block_lines.append(line)

            if line_stripped == "ENDMDL" or line_stripped == "END":
                mol_index += 1
                name = current_name or f"mol_{mol_index:05d}"

                if current_smiles:
                    rows.append({"ligand_id": name, "smiles": current_smiles})
                else:
                    # Try RDKit PDB block conversion as fallback
                    block = "".join(mol_block_lines)
                    mol = Chem.MolFromPDBBlock(block, sanitize=False)
                    if mol is not None:
                        try:
                            Chem.SanitizeMol(mol)
                            smi = Chem.MolToSmiles(mol)
                            if smi:
                                rows.append({"ligand_id": name, "smiles": smi})
                        except Exception:
                            pass

                current_smiles = None
                current_name = None
                mol_block_lines = []

    return pd.DataFrame(rows)


def _from_sdf(path: str) -> pd.DataFrame:
    from rdkit import Chem
    from rdkit.Chem import PandasTools

    df = PandasTools.LoadSDF(path, smilesName="smiles", molColName=None)
    if df is None or df.empty:
        print("[!] No molecules found in SDF file.")
        return pd.DataFrame()

    id_col = None
    for candidate in ["ID", "Name", "PUBCHEM_COMPOUND_CID", "_Name"]:
        if candidate in df.columns:
            id_col = candidate
            break

    if id_col:
        df = df.rename(columns={id_col: "ligand_id"})
    else:
        df["ligand_id"] = [f"mol_{i + 1:05d}" for i in range(len(df))]

    return df[["ligand_id", "smiles"]].dropna()


def _from_smiles_file(path: str) -> pd.DataFrame:
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for idx, line in enumerate(f, start=1):
            parts = line.strip().split(None, 1)
            if not parts:
                continue
            smiles = parts[0]
            name = parts[1] if len(parts) > 1 else f"mol_{idx:05d}"
            rows.append({"ligand_id": name, "smiles": smiles})
    return pd.DataFrame(rows)


def _from_csv(path: str) -> pd.DataFrame:
    # Auto-detect separator: sniff the first line to decide between ; and ,
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        first_line = f.readline()
    sep = ";" if first_line.count(";") > first_line.count(",") else ","

    df = pd.read_csv(
        path,
        sep=sep,
        engine="python",
        on_bad_lines="skip",
        encoding="utf-8",
        encoding_errors="replace",
    )

    print(f"\nDetected separator: '{sep}'")
    print(f"Columns found in your CSV: {list(df.columns)}\n")

    # Already in the right shape
    if {"ligand_id", "smiles"}.issubset(set(df.columns)):
        print("[OK] CSV already has ligand_id and smiles columns.")
        return df[["ligand_id", "smiles"]]

    # Try to auto-detect common alternatives (includes ChEMBL export column names)
    smiles_aliases = ["Smiles", "SMILES", "smiles", "canonical_smiles", "CanonicalSMILES", "smi"]
    id_aliases = [
        "Molecule ChEMBL ID",  # ChEMBL export
        "chembl_id",
        "ID", "id", "Name", "name", "CID", "mol_id", "ligand_id", "zinc_id",
    ]

    smiles_col = next((c for c in smiles_aliases if c in df.columns), None)
    id_col = next((c for c in id_aliases if c in df.columns), None)

    if smiles_col is None:
        print("ERROR: Cannot find a SMILES column. Columns present:", list(df.columns))
        print("       Rename your SMILES column to 'smiles' and rerun.")
        sys.exit(1)

    df = df.rename(columns={smiles_col: "smiles"})

    if id_col is None:
        df["ligand_id"] = [f"mol_{i + 1:05d}" for i in range(len(df))]
    else:
        df = df.rename(columns={id_col: "ligand_id"})

    return df[["ligand_id", "smiles"]]


def convert(input_path: str, output_path: str) -> None:
    ext = os.path.splitext(input_path)[1].lower()

    print(f"[+] Reading: {input_path}")

    if ext == ".pdbqt":
        df = _from_pdbqt(input_path)
    elif ext == ".sdf":
        df = _from_sdf(input_path)
    elif ext in {".smi", ".txt"}:
        df = _from_smiles_file(input_path)
    elif ext == ".csv":
        df = _from_csv(input_path)
    else:
        print(f"[!] Unrecognised extension '{ext}'. Trying CSV parser...")
        df = _from_csv(input_path)

    if df.empty:
        print("ERROR: No molecules could be extracted.")
        sys.exit(1)

    df = df.dropna(subset=["smiles"]).reset_index(drop=True)
    df.to_csv(output_path, index=False)
    print(f"[+] Done. {len(df)} molecules written to: {output_path}")
    print(f"    Columns: {list(df.columns)}")
    print(f"    First row preview: {df.iloc[0].to_dict()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert ligand files to ligands.csv")
    parser.add_argument("--input", required=True, help="Input file (.pdbqt, .sdf, .csv, .smi)")
    parser.add_argument("--output", default="ligands.csv", help="Output CSV path (default: ligands.csv)")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: File not found: {args.input}")
        sys.exit(1)

    convert(args.input, args.output)
