import argparse
import os

import pandas as pd

from clustering import cluster_hits_from_distance
from docking import run_docking
from fingerprints import generate_morgan_fingerprints, tanimoto_similarity_and_distance
from hit_selection import select_cluster_representatives
from ligand_prep import prepare_ligands
from protein_prep import prepare_protein
from visualization import (
    save_chemical_space_pca,
    save_cluster_distribution,
    save_docking_score_distribution,
)


def run_pipeline(
    protein_path: str,
    ligands_csv_path: str,
    output_dir: str,
    top_n: int = 5000,
    max_ligands: int | None = None,
) -> dict:
    os.makedirs(output_dir, exist_ok=True)
    ligand_pdbqt_dir = os.path.join(output_dir, "ligands_pdbqt")
    docking_dir = os.path.join(output_dir, "docking")
    os.makedirs(ligand_pdbqt_dir, exist_ok=True)
    os.makedirs(docking_dir, exist_ok=True)

    docking_results_path = os.path.join(output_dir, "docking_results.csv")
    top_hits_path = os.path.join(output_dir, "top_hits.csv")
    final_hits_path = os.path.join(output_dir, "final_hits.csv")

    docking_plot = os.path.join(output_dir, "docking_score_distribution.png")
    cluster_plot = os.path.join(output_dir, "cluster_distribution.png")
    pca_plot = os.path.join(output_dir, "chemical_space_pca.png")

    print("Loading ligand dataset")
    ligands_df = pd.read_csv(ligands_csv_path)
    print(f"Loaded molecules: {len(ligands_df)}")

    print("Preparing protein")
    receptor_pdbqt = prepare_protein(protein_path, output_dir, receptor_name="protein")

    print("Preparing molecules")
    prepared_df = prepare_ligands(ligands_csv_path, ligand_pdbqt_dir, max_ligands=max_ligands)
    if prepared_df.empty:
        raise RuntimeError("No ligands prepared for docking.")

    print("Running docking")
    docking_df = run_docking(
        receptor_pdbqt=receptor_pdbqt,
        ligands_df=prepared_df,
        output_dir=docking_dir,
        center={"x": 0.0, "y": 0.0, "z": 0.0},
        size={"x": 20.0, "y": 20.0, "z": 20.0},
    )
    docking_df.to_csv(docking_results_path, index=False)

    print("Filtering top molecules")
    top_n = min(top_n, len(docking_df))
    top_hits_df = docking_df.sort_values("docking_score", ascending=True).head(top_n).reset_index(drop=True)
    top_hits_df.to_csv(top_hits_path, index=False)

    print("Generating fingerprints")
    fp_df, fp_matrix = generate_morgan_fingerprints(top_hits_df, radius=2, n_bits=1024)

    print("Calculating Tanimoto similarity and distance")
    _, distance_matrix = tanimoto_similarity_and_distance(fp_matrix)

    print("Running clustering")
    clustered_df = cluster_hits_from_distance(fp_df, distance_matrix, eps=0.7, min_samples=2)

    print("Selecting final hits")
    final_hits_df = select_cluster_representatives(clustered_df)
    final_hits_df.to_csv(final_hits_path, index=False)

    print("Saving results")
    save_docking_score_distribution(docking_df, docking_plot)
    save_cluster_distribution(clustered_df, cluster_plot)
    save_chemical_space_pca(fp_matrix, clustered_df["cluster"], pca_plot)

    return {
        "docking_results": docking_results_path,
        "top_hits": top_hits_path,
        "final_hits": final_hits_path,
        "docking_plot": docking_plot,
        "cluster_plot": cluster_plot,
        "pca_plot": pca_plot,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run virtual screening pipeline")
    parser.add_argument("--protein", default="protein.pdb")
    parser.add_argument("--ligands", default="ligands.csv")
    parser.add_argument("--output-dir", default="data/outputs")
    parser.add_argument("--top-n", type=int, default=5000)
    parser.add_argument("--max-ligands", type=int, default=0)
    args = parser.parse_args()

    max_ligands = args.max_ligands if args.max_ligands > 0 else None
    outputs = run_pipeline(
        protein_path=args.protein,
        ligands_csv_path=args.ligands,
        output_dir=args.output_dir,
        top_n=args.top_n,
        max_ligands=max_ligands,
    )

    print("Pipeline completed")
    for name, path in outputs.items():
        print(f"{name}: {path}")
