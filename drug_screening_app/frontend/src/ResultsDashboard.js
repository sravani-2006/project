import React from "react";

const API_BASE = "http://localhost:8000";

function ClusterBarChart({ data }) {
  if (!data || data.length === 0) {
    return <p className="status">No cluster distribution available yet.</p>;
  }

  const maxSize = Math.max(...data.map((d) => d.size), 1);

  return (
    <div className="chart-wrap">
      {data.map((item) => (
        <div key={item.cluster} className="bar-row">
          <div className="bar-label">Cluster {item.cluster}</div>
          <div className="bar-track">
            <div
              className="bar-fill"
              style={{ width: `${Math.max(6, (item.size / maxSize) * 100)}%` }}
            >
              {item.size}
            </div>
          </div>
        </div>
      ))}
    </div>
  );
}

function ResultsDashboard({ results, pipelineStatus }) {
  return (
    <section className="panel results-panel">
      <h2>Results Dashboard</h2>
      <p className="status-line">Latest status: {pipelineStatus || "Waiting for pipeline"}</p>

      {results?.status !== "completed" ? (
        <p className="status">Run the full pipeline to see results.</p>
      ) : (
        <>
          <div className="stats-grid">
            <div className="stat-card">
              <h3>Total Molecules Screened</h3>
              <p>{results.total_molecules_screened}</p>
            </div>
            <div className="stat-card">
              <h3>Top Molecules Selected</h3>
              <p>{results.top_molecules_selected}</p>
            </div>
            <div className="stat-card">
              <h3>Number of Clusters</h3>
              <p>{results.number_of_clusters}</p>
            </div>
            <div className="stat-card">
              <h3>Final Hit Molecules</h3>
              <p>{results.final_hit_molecules}</p>
            </div>
          </div>

          <h3>Cluster Size Distribution</h3>
          <ClusterBarChart data={results.cluster_distribution} />

          <div className="download-row">
            <a href={`${API_BASE}/download-hits`} className="download-btn" target="_blank" rel="noreferrer">
              Download final_hits.csv
            </a>
            {results.cluster_plot_download ? (
              <a href={`${API_BASE}${results.cluster_plot_download}`} className="download-btn" target="_blank" rel="noreferrer">
                Download cluster plot PNG
              </a>
            ) : null}
          </div>

          <h3>Final Hits Table</h3>
          <div className="table-wrap">
            <table>
              <thead>
                <tr>
                  <th>Ligand ID</th>
                  <th>SMILES</th>
                  <th>Docking Score</th>
                  <th>Cluster</th>
                </tr>
              </thead>
              <tbody>
                {results.final_hits?.map((hit, idx) => (
                  <tr key={`${hit.ligand_id}-${idx}`}>
                    <td>{hit.ligand_id}</td>
                    <td className="smiles-cell">{hit.smiles}</td>
                    <td>{hit.docking_score}</td>
                    <td>{hit.cluster}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </>
      )}
    </section>
  );
}

export default ResultsDashboard;
