import os
import shutil
from datetime import datetime
from typing import Any, Dict, List

import pandas as pd
from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse

from clustering import cluster_hits_from_distance
from docking import run_docking
from fingerprints import generate_morgan_fingerprints, tanimoto_similarity_and_distance
from hit_selection import filter_top_by_docking_score, select_cluster_representatives
from ligand_prep import prepare_ligands
from protein_prep import prepare_protein
from strategy_evaluation import evaluate_strategies
from strategy_selection import (
    select_greedy_diversity_hits,
    select_multiobjective_hits,
    select_score_only_hits,
)
from visualization import (
    save_chemical_space_pca,
    save_cluster_distribution,
    save_docking_score_distribution,
    save_strategy_chemical_space_pca,
    save_strategy_metric_comparison,
)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
UPLOAD_DIR = os.path.join(DATA_DIR, "uploads")
OUTPUT_DIR = os.path.join(DATA_DIR, "outputs")
LIGAND_PDBQT_DIR = os.path.join(DATA_DIR, "ligands_pdbqt")
DOCKING_DIR = os.path.join(DATA_DIR, "docking")

for folder in [DATA_DIR, UPLOAD_DIR, OUTPUT_DIR, LIGAND_PDBQT_DIR, DOCKING_DIR]:
    os.makedirs(folder, exist_ok=True)

app = FastAPI(title="Virtual Screening and Clustering API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

pipeline_state: Dict[str, Any] = {
    "protein_path": None,
    "ligands_path": None,
    "receptor_pdbqt": None,
    "docking_results_path": None,
    "top_hits_path": None,
    "clustered_hits_path": None,
    "final_hits_path": None,
    "strategy1_path": None,
    "strategy2_path": None,
    "strategy3_path": None,
    "strategy4_path": None,
    "strategy_evaluation_path": None,
    "cluster_plot_path": None,
    "docking_plot_path": None,
    "pca_plot_path": None,
    "strategy_score_plot_path": None,
    "strategy_diversity_plot_path": None,
    "strategy_pca_plot_path": None,
    "progress": [],
    "last_updated": None,
}


def log_progress(message: str) -> None:
    print(message)
    pipeline_state["progress"].append(message)
    pipeline_state["last_updated"] = datetime.utcnow().isoformat()


def _safe_remove(path: str) -> None:
    if os.path.exists(path):
        os.remove(path)


def _save_upload_file(upload: UploadFile, destination: str) -> None:
    with open(destination, "wb") as out:
        shutil.copyfileobj(upload.file, out)


def _ensure_input_files() -> None:
    if not pipeline_state["protein_path"] or not os.path.exists(pipeline_state["protein_path"]):
        raise HTTPException(status_code=400, detail="Protein file not uploaded yet.")
    if not pipeline_state["ligands_path"] or not os.path.exists(pipeline_state["ligands_path"]):
        raise HTTPException(status_code=400, detail="Ligand CSV not uploaded yet.")


@app.post("/upload-protein")
async def upload_protein(file: UploadFile = File(...)):
    if not file.filename.lower().endswith(".pdb"):
        raise HTTPException(status_code=400, detail="Protein file must be a .pdb file")

    destination = os.path.join(UPLOAD_DIR, "protein.pdb")
    _save_upload_file(file, destination)

    pipeline_state["protein_path"] = destination
    pipeline_state["progress"] = ["Protein uploaded"]
    pipeline_state["last_updated"] = datetime.utcnow().isoformat()

    return {"message": "Protein uploaded successfully", "path": destination}


def _normalise_ligand_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remap common alternative column names to the required ligand_id / smiles names.
    Handles ChEMBL exports and other common formats.
    """
    col_set = set(df.columns)

    SMILES_ALIASES = ["Smiles", "SMILES", "smiles", "canonical_smiles", "CanonicalSMILES", "smi"]
    ID_ALIASES = [
        "Molecule ChEMBL ID", "chembl_id", "ChEMBL ID",
        "ID", "id", "Name", "name", "CID", "mol_id", "ligand_id", "zinc_id",
    ]

    if "smiles" not in col_set:
        smiles_col = next((c for c in SMILES_ALIASES if c in col_set), None)
        if smiles_col:
            df = df.rename(columns={smiles_col: "smiles"})

    if "ligand_id" not in col_set:
        id_col = next((c for c in ID_ALIASES if c in col_set), None)
        if id_col:
            df = df.rename(columns={id_col: "ligand_id"})
        else:
            df["ligand_id"] = [f"mol_{i + 1:05d}" for i in range(len(df))]

    return df


@app.post("/upload-ligands")
async def upload_ligands(file: UploadFile = File(...)):
    if not file.filename.lower().endswith(".csv"):
        raise HTTPException(status_code=400, detail="Ligand file must be a .csv file")

    destination = os.path.join(UPLOAD_DIR, "ligands.csv")
    _save_upload_file(file, destination)

    try:
        # Auto-detect separator: ChEMBL exports use ; while standard CSVs use ,
        with open(destination, "r", encoding="utf-8", errors="replace") as _f:
            first_line = _f.readline()
        sep = ";" if first_line.count(";") > first_line.count(",") else ","

        ligands_df = pd.read_csv(
            destination,
            sep=sep,
            engine="python",
            on_bad_lines="skip",
            encoding="utf-8",
            encoding_errors="replace",
        )
    except Exception as exc:
        raise HTTPException(
            status_code=400,
            detail=f"Could not parse CSV file: {exc}",
        )

    # Remap ChEMBL / alternative column names to ligand_id and smiles
    ligands_df = _normalise_ligand_df(ligands_df)

    required = {"ligand_id", "smiles"}
    if not required.issubset(set(ligands_df.columns)):
        raise HTTPException(
            status_code=400,
            detail=(
                f"Could not find SMILES or ID columns. "
                f"Columns detected: {list(ligands_df.columns)}. "
                f"Expected a column named one of: Smiles, SMILES, canonical_smiles "
                f"and an ID column like 'Molecule ChEMBL ID', ID, Name, CID."
            ),
        )

    # Drop rows with missing SMILES and save normalised file in place
    ligands_df = ligands_df.dropna(subset=["smiles"])
    ligands_df[["ligand_id", "smiles"]].to_csv(destination, index=False)

    pipeline_state["ligands_path"] = destination
    pipeline_state["progress"].append("Ligands uploaded")
    pipeline_state["last_updated"] = datetime.utcnow().isoformat()

    return {
        "message": "Ligands uploaded successfully",
        "path": destination,
        "total_rows": int(len(ligands_df)),
    }


@app.post("/run-docking")
async def run_docking_endpoint(max_ligands: int = 50):
    _ensure_input_files()

    log_progress("Preparing protein")
    receptor_pdbqt = prepare_protein(pipeline_state["protein_path"], OUTPUT_DIR, receptor_name="protein")
    pipeline_state["receptor_pdbqt"] = receptor_pdbqt

    log_progress("Preparing ligands")
    # Allow UI-driven run size while keeping safe defaults for demo responsiveness.
    run_cap = max(0, int(max_ligands))
    run_cap = min(run_cap, 5000)
    max_ligands_arg = run_cap if run_cap > 0 else None
    prepared_ligands_df = prepare_ligands(
        pipeline_state["ligands_path"],
        LIGAND_PDBQT_DIR,
        max_ligands=max_ligands_arg,
    )
    if prepared_ligands_df.empty:
        raise HTTPException(status_code=400, detail="No ligands were successfully prepared.")

    log_progress("Running docking")
    # Default search box center/size can be updated once the user knows target pocket coordinates.
    docking_df = run_docking(
        receptor_pdbqt=receptor_pdbqt,
        ligands_df=prepared_ligands_df,
        output_dir=DOCKING_DIR,
        center={"x": 0.0, "y": 0.0, "z": 0.0},
        size={"x": 20.0, "y": 20.0, "z": 20.0},
    )

    docking_results_path = os.path.join(OUTPUT_DIR, "docking_results.csv")
    docking_df.to_csv(docking_results_path, index=False)
    pipeline_state["docking_results_path"] = docking_results_path

    docking_plot_path = save_docking_score_distribution(
        docking_df, os.path.join(OUTPUT_DIR, "docking_score_distribution.png")
    )
    pipeline_state["docking_plot_path"] = docking_plot_path

    return {
        "message": "Docking completed",
        "docked_molecules": int(len(docking_df)),
        "output": docking_results_path,
    }


@app.post("/filter-hits")
async def filter_hits_endpoint(top_n: int = 5000):
    path = pipeline_state["docking_results_path"]
    if not path or not os.path.exists(path):
        raise HTTPException(status_code=400, detail="Run docking first.")

    log_progress("Filtering hits")
    docking_df = pd.read_csv(path)
    if docking_df.empty:
        raise HTTPException(status_code=400, detail="Docking results are empty.")

    top_hits_df = filter_top_by_docking_score(docking_df, top_n=top_n)
    top_hits_path = os.path.join(OUTPUT_DIR, "top_hits.csv")
    top_hits_df.to_csv(top_hits_path, index=False)

    pipeline_state["top_hits_path"] = top_hits_path

    return {
        "message": "Top hits generated",
        "selected": int(len(top_hits_df)),
        "output": top_hits_path,
    }


@app.post("/run-clustering")
async def run_clustering_endpoint(eps: float = 0.7, min_samples: int = 2):
    top_hits_path = pipeline_state["top_hits_path"]
    if not top_hits_path or not os.path.exists(top_hits_path):
        raise HTTPException(status_code=400, detail="Run /filter-hits first.")

    top_hits_df = pd.read_csv(top_hits_path)
    if top_hits_df.empty:
        raise HTTPException(status_code=400, detail="Top hits are empty.")

    log_progress("Generating fingerprints")
    filtered_df, fp_matrix = generate_morgan_fingerprints(top_hits_df, radius=2, n_bits=1024)
    if len(filtered_df) == 0:
        raise HTTPException(status_code=400, detail="No valid molecules for fingerprint generation.")

    log_progress("Calculating similarity")
    _, distance_matrix = tanimoto_similarity_and_distance(fp_matrix)

    log_progress("Clustering molecules")
    clustered_df = cluster_hits_from_distance(filtered_df, distance_matrix, eps=eps, min_samples=min_samples)
    clustered_hits_path = os.path.join(OUTPUT_DIR, "clustered_hits.csv")
    clustered_df.to_csv(clustered_hits_path, index=False)
    pipeline_state["clustered_hits_path"] = clustered_hits_path

    log_progress("Selecting diverse hits")
    final_hits_df = select_cluster_representatives(clustered_df)
    final_hits_path = os.path.join(OUTPUT_DIR, "final_hits.csv")
    final_hits_df.to_csv(final_hits_path, index=False)
    pipeline_state["final_hits_path"] = final_hits_path

    cluster_plot_path = save_cluster_distribution(clustered_df, os.path.join(OUTPUT_DIR, "cluster_distribution.png"))
    pipeline_state["cluster_plot_path"] = cluster_plot_path

    pca_plot_path = save_chemical_space_pca(
        fp_matrix,
        clustered_df["cluster"],
        os.path.join(OUTPUT_DIR, "chemical_space_pca.png"),
    )
    pipeline_state["pca_plot_path"] = pca_plot_path

    n_clusters = int(clustered_df[clustered_df["cluster"] != -1]["cluster"].nunique())

    return {
        "message": "Clustering and representative hit selection completed",
        "clusters": n_clusters,
        "final_hits": int(len(final_hits_df)),
        "clustered_output": clustered_hits_path,
        "final_hits_output": final_hits_path,
    }


@app.post("/run-strategy-analysis")
async def run_strategy_analysis_endpoint(k: int = 20):
    docking_path = pipeline_state["docking_results_path"]
    if not docking_path or not os.path.exists(docking_path):
        raise HTTPException(status_code=400, detail="Run docking first.")

    k = max(1, int(k))

    # Ensure DBSCAN strategy artifacts exist for fair comparison.
    if not pipeline_state["top_hits_path"] or not os.path.exists(pipeline_state["top_hits_path"]):
        await filter_hits_endpoint()
    if not pipeline_state["final_hits_path"] or not os.path.exists(pipeline_state["final_hits_path"]):
        await run_clustering_endpoint()

    log_progress("Running Strategy 1 (Score Only)")
    docking_df = pd.read_csv(docking_path)
    strategy1_df = select_score_only_hits(docking_df, k=k)
    strategy1_path = os.path.join(OUTPUT_DIR, "strategy1_score_only_hits.csv")
    strategy1_df.to_csv(strategy1_path, index=False)
    pipeline_state["strategy1_path"] = strategy1_path

    log_progress("Preparing Strategy 2 (DBSCAN Selection)")
    strategy2_source = pd.read_csv(pipeline_state["final_hits_path"])
    strategy2_df = strategy2_source[["ligand_id", "smiles", "docking_score", "cluster"]].copy()
    strategy2_df = strategy2_df.sort_values("docking_score", ascending=True).head(k).reset_index(drop=True)
    strategy2_path = os.path.join(OUTPUT_DIR, "strategy2_dbscan_hits.csv")
    strategy2_df.to_csv(strategy2_path, index=False)
    pipeline_state["strategy2_path"] = strategy2_path

    log_progress("Running Strategy 3 (Greedy Diversity)")
    strategy3_df = select_greedy_diversity_hits(docking_df, k=k)
    strategy3_path = os.path.join(OUTPUT_DIR, "strategy3_greedy_diversity_hits.csv")
    strategy3_df.to_csv(strategy3_path, index=False)
    pipeline_state["strategy3_path"] = strategy3_path

    log_progress("Running Strategy 4 (Multi-Objective)")
    strategy4_df = select_multiobjective_hits(docking_df, k=k)
    strategy4_path = os.path.join(OUTPUT_DIR, "strategy4_multiobjective_hits.csv")
    strategy4_df.to_csv(strategy4_path, index=False)
    pipeline_state["strategy4_path"] = strategy4_path

    log_progress("Evaluating strategies")
    strategy_to_hits = {
        "Strategy 1: Score Only": strategy1_df,
        "Strategy 2: DBSCAN": strategy2_df,
        "Strategy 3: Greedy Diversity": strategy3_df,
        "Strategy 4: Multi-Objective": strategy4_df,
    }
    evaluation_df = evaluate_strategies(strategy_to_hits)
    evaluation_path = os.path.join(OUTPUT_DIR, "strategy_evaluation.csv")
    evaluation_df.to_csv(evaluation_path, index=False)
    pipeline_state["strategy_evaluation_path"] = evaluation_path

    score_plot_path = save_strategy_metric_comparison(
        evaluation_df,
        metric_column="average_docking_score",
        title="Average Docking Score by Strategy",
        y_label="Average Docking Score",
        output_path=os.path.join(OUTPUT_DIR, "strategy_score_comparison.png"),
    )
    diversity_plot_path = save_strategy_metric_comparison(
        evaluation_df,
        metric_column="average_similarity",
        title="Average Similarity by Strategy",
        y_label="Average Pairwise Similarity",
        output_path=os.path.join(OUTPUT_DIR, "strategy_diversity_comparison.png"),
    )
    strategy_pca_plot_path = save_strategy_chemical_space_pca(
        strategy_to_hits,
        output_path=os.path.join(OUTPUT_DIR, "chemical_space_pca.png"),
    )

    pipeline_state["strategy_score_plot_path"] = score_plot_path
    pipeline_state["strategy_diversity_plot_path"] = diversity_plot_path
    pipeline_state["strategy_pca_plot_path"] = strategy_pca_plot_path

    return {
        "message": "Strategy analysis completed",
        "selected_per_strategy": k,
        "strategy1_output": strategy1_path,
        "strategy2_output": strategy2_path,
        "strategy3_output": strategy3_path,
        "strategy4_output": strategy4_path,
        "evaluation_output": evaluation_path,
    }


@app.post("/run-full-pipeline")
async def run_full_pipeline():
    await run_docking_endpoint()
    await filter_hits_endpoint()
    await run_clustering_endpoint()
    await run_strategy_analysis_endpoint()

    return {
        "message": "Full pipeline completed",
        "results_endpoint": "/results",
        "download_endpoint": "/download-hits",
    }


@app.get("/results")
async def get_results():
    docking_path = pipeline_state["docking_results_path"]
    top_hits_path = pipeline_state["top_hits_path"]
    clustered_path = pipeline_state["clustered_hits_path"]
    final_path = pipeline_state["final_hits_path"]
    strategy1_path = pipeline_state["strategy1_path"]
    strategy2_path = pipeline_state["strategy2_path"]
    strategy3_path = pipeline_state["strategy3_path"]
    strategy4_path = pipeline_state["strategy4_path"]
    evaluation_path = pipeline_state["strategy_evaluation_path"]

    if not final_path or not os.path.exists(final_path):
        return {
            "status": "pending",
            "progress": pipeline_state["progress"],
            "message": "Pipeline not completed yet",
        }

    docking_df = pd.read_csv(docking_path) if docking_path and os.path.exists(docking_path) else pd.DataFrame()
    top_hits_df = pd.read_csv(top_hits_path) if top_hits_path and os.path.exists(top_hits_path) else pd.DataFrame()
    clustered_df = pd.read_csv(clustered_path) if clustered_path and os.path.exists(clustered_path) else pd.DataFrame()
    final_df = pd.read_csv(final_path)
    strategy1_df = pd.read_csv(strategy1_path) if strategy1_path and os.path.exists(strategy1_path) else pd.DataFrame()
    strategy2_df = pd.read_csv(strategy2_path) if strategy2_path and os.path.exists(strategy2_path) else pd.DataFrame()
    strategy3_df = pd.read_csv(strategy3_path) if strategy3_path and os.path.exists(strategy3_path) else pd.DataFrame()
    strategy4_df = pd.read_csv(strategy4_path) if strategy4_path and os.path.exists(strategy4_path) else pd.DataFrame()
    evaluation_df = pd.read_csv(evaluation_path) if evaluation_path and os.path.exists(evaluation_path) else pd.DataFrame()
    if not evaluation_df.empty:
        evaluation_df = evaluation_df.where(pd.notna(evaluation_df), None)

    cluster_distribution: List[Dict[str, Any]] = []
    if not clustered_df.empty:
        counts = clustered_df[clustered_df["cluster"] != -1]["cluster"].value_counts().sort_index()
        cluster_distribution = [
            {"cluster": int(cluster_id), "size": int(size)} for cluster_id, size in counts.items()
        ]

    return {
        "status": "completed",
        "progress": pipeline_state["progress"],
        "total_molecules_screened": int(len(docking_df)),
        "top_molecules_selected": int(len(top_hits_df)),
        "number_of_clusters": int(clustered_df[clustered_df["cluster"] != -1]["cluster"].nunique()) if not clustered_df.empty else 0,
        "final_hit_molecules": int(len(final_df)),
        "cluster_distribution": cluster_distribution,
        "docking_results": docking_df.to_dict(orient="records"),
        "top_hits": top_hits_df.to_dict(orient="records"),
        "clustered_hits": clustered_df.to_dict(orient="records") if not clustered_df.empty else [],
        "final_hits": final_df.to_dict(orient="records"),
        "strategy1_hits": strategy1_df.to_dict(orient="records") if not strategy1_df.empty else [],
        "strategy2_hits": strategy2_df.to_dict(orient="records") if not strategy2_df.empty else [],
        "strategy3_hits": strategy3_df.to_dict(orient="records") if not strategy3_df.empty else [],
        "strategy4_hits": strategy4_df.to_dict(orient="records") if not strategy4_df.empty else [],
        "strategy_evaluation": evaluation_df.to_dict(orient="records") if not evaluation_df.empty else [],
        "strategy_comparison_table": evaluation_df.to_dict(orient="records") if not evaluation_df.empty else [],
        "docking_plot_download": "/download-plot/docking" if pipeline_state["docking_plot_path"] else None,
        "cluster_plot_download": "/download-plot/cluster" if pipeline_state["cluster_plot_path"] else None,
        "pca_plot_download": "/download-plot/pca" if pipeline_state["pca_plot_path"] else None,
        "strategy_score_plot_download": "/download-plot/strategy-score" if pipeline_state["strategy_score_plot_path"] else None,
        "strategy_diversity_plot_download": "/download-plot/strategy-diversity" if pipeline_state["strategy_diversity_plot_path"] else None,
        "strategy_pca_plot_download": "/download-plot/strategy-pca" if pipeline_state["strategy_pca_plot_path"] else None,
        "strategy1_download": "/download-strategy/score-only" if pipeline_state["strategy1_path"] else None,
        "strategy2_download": "/download-strategy/dbscan" if pipeline_state["strategy2_path"] else None,
        "strategy3_download": "/download-strategy/greedy-diversity" if pipeline_state["strategy3_path"] else None,
        "strategy4_download": "/download-strategy/multi-objective" if pipeline_state["strategy4_path"] else None,
        "strategy_evaluation_download": "/download-evaluation" if pipeline_state["strategy_evaluation_path"] else None,
    }


@app.get("/download-hits")
async def download_hits():
    file_path = pipeline_state["final_hits_path"]
    if not file_path or not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="No final hits available yet.")
    return FileResponse(file_path, media_type="text/csv", filename="final_hits.csv")


@app.get("/download-strategy/{strategy_id}")
async def download_strategy(strategy_id: str):
    mapping = {
        "score-only": (pipeline_state["strategy1_path"], "strategy1_score_only_hits.csv"),
        "dbscan": (pipeline_state["strategy2_path"], "strategy2_dbscan_hits.csv"),
        "greedy-diversity": (pipeline_state["strategy3_path"], "strategy3_greedy_diversity_hits.csv"),
        "multi-objective": (pipeline_state["strategy4_path"], "strategy4_multiobjective_hits.csv"),
    }
    if strategy_id not in mapping:
        raise HTTPException(status_code=404, detail="Unknown strategy")

    file_path, filename = mapping[strategy_id]
    if not file_path or not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail=f"No output for strategy '{strategy_id}' yet.")
    return FileResponse(file_path, media_type="text/csv", filename=filename)


@app.get("/download-evaluation")
async def download_evaluation():
    file_path = pipeline_state["strategy_evaluation_path"]
    if not file_path or not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="No strategy evaluation available yet.")
    return FileResponse(file_path, media_type="text/csv", filename="strategy_evaluation.csv")


@app.get("/download-plot/{plot_kind}")
async def download_plot(plot_kind: str):
    mapping = {
        "docking": (pipeline_state["docking_plot_path"], "docking_score_distribution.png"),
        "cluster": (pipeline_state["cluster_plot_path"], "cluster_distribution.png"),
        "pca": (pipeline_state["pca_plot_path"], "chemical_space_pca.png"),
        "strategy-score": (pipeline_state["strategy_score_plot_path"], "strategy_score_comparison.png"),
        "strategy-diversity": (pipeline_state["strategy_diversity_plot_path"], "strategy_diversity_comparison.png"),
        "strategy-pca": (pipeline_state["strategy_pca_plot_path"], "chemical_space_pca.png"),
    }
    if plot_kind not in mapping:
        raise HTTPException(status_code=404, detail="Unknown plot type")

    file_path, filename = mapping[plot_kind]
    if not file_path or not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail=f"No {plot_kind} plot available yet.")
    return FileResponse(file_path, media_type="image/png", filename=filename)
