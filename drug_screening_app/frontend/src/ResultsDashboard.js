import React, { useMemo, useState } from "react";

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

function ResultsTable({ rows, columns }) {
  return (
    <div className="table-wrap">
      <table>
        <thead>
          <tr>
            {columns.map((col) => (
              <th key={col.key}>{col.label}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {rows?.slice(0, 500).map((row, idx) => (
            <tr key={`${idx}-${row.ligand_id || "row"}`}>
              {columns.map((col) => (
                <td key={col.key} className={col.key === "smiles" ? "smiles-cell" : ""}>
                  {String(row[col.key] ?? "")}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function ResultsDashboard({ results, pipelineStatus }) {
  const [activeTab, setActiveTab] = useState("final");

  const tabs = useMemo(
    () => [
      { id: "docking", label: "Docking Results" },
      { id: "top", label: "Top Hits" },
      { id: "clusters", label: "Clusters" },
      { id: "final", label: "Final Candidates" },
      { id: "viz", label: "Visualizations" },
    ],
    []
  );

  const columns = [
    { key: "ligand_id", label: "Ligand ID" },
    { key: "smiles", label: "SMILES" },
    { key: "docking_score", label: "Docking Score" },
    { key: "cluster", label: "Cluster" },
  ];

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
            {results.docking_plot_download ? (
              <a href={`${API_BASE}${results.docking_plot_download}`} className="download-btn" target="_blank" rel="noreferrer">
                Download docking score PNG
              </a>
            ) : null}
            {results.cluster_plot_download ? (
              <a href={`${API_BASE}${results.cluster_plot_download}`} className="download-btn" target="_blank" rel="noreferrer">
                Download cluster plot PNG
              </a>
            ) : null}
            {results.pca_plot_download ? (
              <a href={`${API_BASE}${results.pca_plot_download}`} className="download-btn" target="_blank" rel="noreferrer">
                Download chemical space PCA PNG
              </a>
            ) : null}
          </div>

          <div className="tab-row">
            {tabs.map((tab) => (
              <button
                key={tab.id}
                className={`tab-btn ${activeTab === tab.id ? "tab-active" : ""}`}
                onClick={() => setActiveTab(tab.id)}
              >
                {tab.label}
              </button>
            ))}
          </div>

          {activeTab === "docking" ? <ResultsTable rows={results.docking_results || []} columns={columns} /> : null}
          {activeTab === "top" ? <ResultsTable rows={results.top_hits || []} columns={columns} /> : null}
          {activeTab === "clusters" ? <ResultsTable rows={results.clustered_hits || []} columns={columns} /> : null}
          {activeTab === "final" ? <ResultsTable rows={results.final_hits || []} columns={columns} /> : null}

          {activeTab === "viz" ? (
            <div className="viz-grid">
              {results.docking_plot_download ? (
                <div className="viz-card">
                  <h4>Docking Score Histogram</h4>
                  <img src={`${API_BASE}${results.docking_plot_download}`} alt="Docking score distribution" />
                </div>
              ) : null}
              {results.cluster_plot_download ? (
                <div className="viz-card">
                  <h4>Cluster Size Distribution</h4>
                  <img src={`${API_BASE}${results.cluster_plot_download}`} alt="Cluster size distribution" />
                </div>
              ) : null}
              {results.pca_plot_download ? (
                <div className="viz-card">
                  <h4>Chemical Space PCA</h4>
                  <img src={`${API_BASE}${results.pca_plot_download}`} alt="Chemical space PCA" />
                </div>
              ) : null}
            </div>
          ) : null}
        </>
      )}
    </section>
  );
}

export default ResultsDashboard;
