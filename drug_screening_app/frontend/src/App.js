import React, { useCallback, useEffect, useState } from "react";
import UploadProtein from "./UploadProtein";
import UploadLigands from "./UploadLigands";
import RunPipeline from "./RunPipeline";
import ResultsDashboard from "./ResultsDashboard";
import "./styles.css";

const API_BASE = "http://localhost:8000";

function App() {
  const [pipelineStatus, setPipelineStatus] = useState("");
  const [results, setResults] = useState(null);

  const fetchResults = useCallback(async () => {
    try {
      const response = await fetch(`${API_BASE}/results`);
      const data = await response.json();
      setResults(data);

      if (data.progress && data.progress.length > 0) {
        setPipelineStatus(data.progress[data.progress.length - 1]);
      }
    } catch (err) {
      setPipelineStatus(`Failed to fetch results: ${err.message}`);
    }
  }, []);

  useEffect(() => {
    fetchResults();
  }, [fetchResults]);

  return (
    <div className="app-shell">
      <header className="hero">
        <div className="hero-content">
          <h1>Virtual Screening and Molecule Clustering</h1>
          <p>
            Upload a protein and ligand library, run docking, cluster top hits, and select diverse
            candidate molecules for drug discovery.
          </p>
        </div>
      </header>

      <main className="content-grid">
        <UploadProtein onUploaded={() => setPipelineStatus("Protein uploaded")} />
        <UploadLigands onUploaded={() => setPipelineStatus("Ligands uploaded")} />
        <RunPipeline onPipelineDone={fetchResults} onStatus={setPipelineStatus} />
        <ResultsDashboard results={results} pipelineStatus={pipelineStatus} />
      </main>
    </div>
  );
}

export default App;
