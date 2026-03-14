import React, { useState } from "react";

const API_BASE = "http://localhost:8000";

function RunPipeline({ onPipelineDone, onStatus }) {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");

  const runVirtualScreening = async () => {
    try {
      setLoading(true);
      setError("");

      onStatus("Preparing protein");
      let response = await fetch(`${API_BASE}/run-docking`, { method: "POST" });
      if (!response.ok) throw new Error((await response.json()).detail || "Docking stage failed");

      onStatus("Filtering hits");
      response = await fetch(`${API_BASE}/filter-hits`, { method: "POST" });
      if (!response.ok) throw new Error((await response.json()).detail || "Filtering stage failed");

      onStatus("Generating fingerprints");
      onStatus("Clustering molecules");
      onStatus("Selecting diverse hits");
      response = await fetch(`${API_BASE}/run-clustering`, { method: "POST" });
      if (!response.ok) throw new Error((await response.json()).detail || "Clustering stage failed");

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
      <button onClick={runVirtualScreening} disabled={loading} className="run-button">
        {loading ? "Running..." : "Run Virtual Screening"}
      </button>

      {error ? <p className="error">{error}</p> : null}
    </section>
  );
}

export default RunPipeline;
