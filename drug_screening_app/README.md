# Drug Screening App

A full-stack virtual screening workflow for:

1. Protein preparation
2. Ligand preparation from SMILES
3. AutoDock Vina docking
4. Top-hit filtering
5. Morgan fingerprint generation
6. DBSCAN clustering
7. Diverse representative hit selection

## Project Structure

```text
drug_screening_app/
  backend/
    main.py
    protein_prep.py
    ligand_prep.py
    docking.py
    fingerprints.py
    clustering.py
    hit_selection.py
    requirements.txt
  frontend/
    package.json
    public/
      index.html
    src/
      App.js
      UploadProtein.js
      UploadLigands.js
      RunPipeline.js
      ResultsDashboard.js
      styles.css
      index.js
```

## Backend Setup (FastAPI)

From `drug_screening_app/backend`:

```bash
pip install -r requirements.txt
uvicorn main:app --reload
```

Backend URL: `http://localhost:8000`

## Frontend Setup (React)

From `drug_screening_app/frontend`:

```bash
npm install
npm start
```

Frontend URL: `http://localhost:3000`

## Required Input Files

1. Protein file: `protein.pdb`
2. Ligands file: `ligands.csv` with columns:
   - `ligand_id`
   - `smiles`

## API Endpoints

1. `POST /upload-protein`
2. `POST /upload-ligands`
3. `POST /run-docking`
4. `POST /filter-hits`
5. `POST /run-clustering`
6. `GET /results`
7. `GET /download-hits`

Convenience endpoint:

1. `POST /run-full-pipeline`

## Output Files

Written to `backend/data/outputs/`:

1. `docking_results.csv`
2. `top_hits.csv`
3. `clustered_hits.csv`
4. `final_hits.csv`
5. `cluster_distribution.png`

## Notes on Docking Tools

For real docking, install:

1. AutoDock Vina (`vina` executable)
2. Open Babel (`obabel` executable)

If Vina is unavailable, the app uses deterministic mock docking scores so the complete pipeline and UI remain runnable for learning and testing.
