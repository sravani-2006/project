import os
import shutil
import subprocess
from typing import Optional


def _remove_waters_from_pdb(input_path: str, output_path: str) -> None:
    """Remove water molecules (HOH/WAT) from a PDB file."""
    with open(input_path, "r", encoding="utf-8") as src, open(output_path, "w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith(("ATOM", "HETATM")):
                residue_name = line[17:20].strip().upper()
                if residue_name in {"HOH", "WAT"}:
                    continue
            dst.write(line)


def _run_command(command: list[str]) -> bool:
    try:
        completed = subprocess.run(command, check=True, capture_output=True, text=True)
        if completed.stdout.strip():
            print(completed.stdout)
        if completed.stderr.strip():
            print(completed.stderr)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"Command failed: {' '.join(command)}")
        print(str(exc))
        return False


def prepare_protein(protein_pdb_path: str, output_dir: str, receptor_name: Optional[str] = "protein") -> str:
    """
    Prepare a receptor for docking.

    Steps:
    1) Remove waters
    2) Add hydrogens and convert to PDBQT using common external tools
    """
    os.makedirs(output_dir, exist_ok=True)

    cleaned_pdb = os.path.join(output_dir, f"{receptor_name}_cleaned.pdb")
    receptor_pdbqt = os.path.join(output_dir, f"{receptor_name}.pdbqt")

    print("[Protein Prep] Removing water molecules...")
    _remove_waters_from_pdb(protein_pdb_path, cleaned_pdb)

    print("[Protein Prep] Adding hydrogens and converting to PDBQT...")

    # Try Open Babel first.
    obabel_cmd = ["obabel", cleaned_pdb, "-O", receptor_pdbqt, "-h"]
    if _run_command(obabel_cmd):
        print(f"[Protein Prep] Receptor prepared: {receptor_pdbqt}")
        return receptor_pdbqt

    # Fallback to AutoDockTools script if available.
    adt_cmd = ["prepare_receptor4.py", "-r", cleaned_pdb, "-o", receptor_pdbqt, "-A", "hydrogens"]
    if _run_command(adt_cmd):
        print(f"[Protein Prep] Receptor prepared: {receptor_pdbqt}")
        return receptor_pdbqt

    # Keep cleaned file for troubleshooting but fail clearly.
    if os.path.exists(receptor_pdbqt):
        os.remove(receptor_pdbqt)
    raise RuntimeError(
        "Could not convert protein to PDBQT. Install Open Babel (`obabel`) or AutoDockTools (`prepare_receptor4.py`)."
    )
