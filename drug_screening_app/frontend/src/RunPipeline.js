import React, { useState } from "react";

const API_BASE = "http://127.0.0.1:8000";

function RunPipeline({ onPipelineDone, onStatus }) {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");
  const [maxLigands, setMaxLigands] = useState("50");

  const runVirtualScreening = async () => {
    try {
      setLoading(true);
      setError("");

      onStatus("Preparing protein");
      let response = await fetch(`${API_BASE}/run-docking?max_ligands=${maxLigands}`, { method: "POST" });
      if (!response.ok) throw new Error((await response.json()).detail || "Docking stage failed");

      onStatus("Filtering hits");
      response = await fetch(`${API_BASE}/filter-hits`, { method: "POST" });
      if (!response.ok) throw new Error((await response.json()).detail || "Filtering stage failed");

      onStatus("Generating fingerprints");
      onStatus("Clustering molecules");
      onStatus("Selecting diverse hits");
      response = await fetch(`${API_BASE}/run-clustering`, { method: "POST" });
      if (!response.ok) throw new Error((await response.json()).detail || "Clustering stage failed");

      onStatus("Running hit-selection strategies");
      response = await fetch(`${API_BASE}/run-strategy-analysis?k=20`, { method: "POST" });
      if (!response.ok) throw new Error((await response.json()).detail || "Strategy analysis failed");

      onStatus("Pipeline complete");
      onPipelineDone();
    } catch (err) {
      setError(err.message);
      onStatus(`Pipeline error: ${err.message}`);
    } finally {
      setLoading(false);
    }
  };

  return (
    <section className="panel">
      <h2>Run Pipeline</h2>

      <label className="run-size-label" htmlFor="run-size-select">
        Molecules To Process
      </label>
      <select
        id="run-size-select"
        className="run-size-select"
        value={maxLigands}
        onChange={(e) => setMaxLigands(e.target.value)}
        disabled={loading}
      >
        <option value="20">20 (fast)</option>
        <option value="50">50 (recommended)</option>
        <option value="100">100 (more clusters)</option>
        <option value="0">Full uploaded dataset</option>
      </select>

      <button onClick={runVirtualScreening} disabled={loading} className="run-button">
        {loading ? "Running..." : "Run Virtual Screening"}
      </button>

      {error ? <p className="error">{error}</p> : null}
    </section>
  );
}

export default RunPipeline;
