import argparse
import os

import pandas as pd

from clustering import cluster_hits
from docking import run_docking
from fingerprints import generate_morgan_fingerprints
from hit_selection import filter_top_by_docking_score, select_cluster_representatives
from ligand_prep import prepare_ligands
from protein_prep import prepare_protein


def run_demo_pipeline(
    protein_path: str,
    ligands_path: str,
    work_dir: str,
    demo_size: int = 100,
) -> None:
    os.makedirs(work_dir, exist_ok=True)

    ligands_pdbqt_dir = os.path.join(work_dir, "ligands_pdbqt_demo")
    docking_dir = os.path.join(work_dir, "docking_demo")
    os.makedirs(ligands_pdbqt_dir, exist_ok=True)
    os.makedirs(docking_dir, exist_ok=True)

    # Required output files from the user request.
    demo_dataset_path = os.path.join(work_dir, "ligands_demo_100.csv")
    docking_results_path = os.path.join(work_dir, "docking_results_demo.csv")
    top_hits_path = os.path.join(work_dir, "top_hits_demo.csv")
    final_hits_path = os.path.join(work_dir, "final_hits_demo.csv")

    print("Loading ligand dataset")
    ligands_df = pd.read_csv(ligands_path)
    required_cols = {"ligand_id", "smiles"}
    missing = required_cols - set(ligands_df.columns)
    if missing:
        raise ValueError(f"Ligands CSV is missing required columns: {sorted(missing)}")

    print("Selecting 100 molecules for demo")
    n = min(demo_size, len(ligands_df))
    demo_df = ligands_df.sample(n=n, random_state=42).reset_index(drop=True)
    demo_df.to_csv(demo_dataset_path, index=False)
    print(f"Saved demo dataset: {demo_dataset_path} (rows={len(demo_df)})")

    print("Preparing protein")
    receptor_pdbqt = prepare_protein(protein_path, work_dir, receptor_name="protein_demo")

    print("Preparing ligands")
    prepared_ligands_df = prepare_ligands(demo_dataset_path, ligands_pdbqt_dir)
    if prepared_ligands_df.empty:
        raise RuntimeError("No ligands were prepared successfully for demo run.")

    print("Running docking")
    docking_df = run_docking(
        receptor_pdbqt=receptor_pdbqt,
        ligands_df=prepared_ligands_df,
        output_dir=docking_dir,
        center={"x": 0.0, "y": 0.0, "z": 0.0},
        size={"x": 20.0, "y": 20.0, "z": 20.0},
        vina_executable="vina",
        allow_mock_if_vina_missing=True,
    )
    docking_df.to_csv(docking_results_path, index=False)
    print(f"Saved docking results: {docking_results_path}")

    print("Filtering top molecules")
    top_hits_df = filter_top_by_docking_score(docking_df, top_n=30)
    top_hits_df.to_csv(top_hits_path, index=False)
    print(f"Saved top hits: {top_hits_path}")

    print("Generating fingerprints")
    fp_input_df, fp_matrix = generate_morgan_fingerprints(top_hits_df, radius=2, n_bits=1024)
    if len(fp_input_df) == 0:
        raise RuntimeError("No valid molecules after fingerprint generation.")

    print("Running clustering")
    clustered_df = cluster_hits(fp_input_df, fp_matrix, eps=0.7, min_samples=2)

    print("Selecting representative molecules")
    final_hits_df = select_cluster_representatives(clustered_df)
    print(f"Number of final hit molecules: {len(final_hits_df)}")
    final_hits_df.to_csv(final_hits_path, index=False)
    print(f"Saved final hits: {final_hits_path}")

    print("Demo pipeline completed")
    print("Generated files:")
    print(f" - {docking_results_path}")
    print(f" - {top_hits_path}")
    print(f" - {final_hits_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a 100-ligand demo virtual screening pipeline")
    parser.add_argument("--protein", default="protein.pdb", help="Path to protein PDB file")
    parser.add_argument("--ligands", default="ligands.csv", help="Path to ligand CSV file")
    parser.add_argument(
        "--work-dir",
        default=os.path.join("data", "demo"),
        help="Directory for demo outputs and intermediate files",
    )
    args = parser.parse_args()

    run_demo_pipeline(
        protein_path=args.protein,
        ligands_path=args.ligands,
        work_dir=args.work_dir,
        demo_size=100,
    )
