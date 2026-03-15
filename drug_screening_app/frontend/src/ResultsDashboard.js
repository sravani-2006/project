import React, { useMemo, useState } from "react";

const API_BASE = "http://127.0.0.1:8000";

function getBindingStrength(score) {
  const value = Number(score);
  if (!Number.isFinite(value)) return "";
  if (value <= -10) return "very strong";
  if (value <= -9) return "strong";
  if (value <= -8) return "good";
  if (value <= -7) return "moderate";
  return "weak";
}

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
              {columns.map((col) => {
                const cellValue = col.render ? col.render(row) : String(row[col.key] ?? "");
                const cellClass = col.className || (col.key === "smiles" ? "smiles-cell" : "");
                return (
                  <td key={col.key} className={cellClass}>
                    {cellValue}
                  </td>
                );
              })}
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
      { id: "score_based", label: "Score-Based Selection" },
      { id: "dbscan_combined", label: "DBSCAN and Clusters" },
      { id: "final", label: "Final Candidates" },
      { id: "greedy", label: "Greedy Diversity" },
      { id: "multi", label: "Multi-Objective" },
      { id: "evaluation", label: "Evaluation" },
      { id: "viz", label: "Visualizations" },
    ],
    []
  );

  const columns = [
    { key: "ligand_id", label: "Ligand ID" },
    { key: "smiles", label: "SMILES" },
    { key: "docking_score", label: "Docking Score" },
    {
      key: "binding_strength",
      label: "Binding Strength",
      className: "strength-cell",
      render: (row) => {
        const strength = getBindingStrength(row.docking_score);
        return strength ? <span className={`strength-chip strength-${strength.replace(" ", "-")}`}>{strength}</span> : "";
      },
    },
    { key: "cluster", label: "Cluster" },
  ];

  const strategyColumns = [
    { key: "ligand_id", label: "Ligand ID" },
    { key: "smiles", label: "SMILES" },
    { key: "docking_score", label: "Docking Score" },
  ];

  const strategy2Columns = [
    { key: "ligand_id", label: "Ligand ID" },
    { key: "smiles", label: "SMILES" },
    { key: "docking_score", label: "Docking Score" },
    { key: "cluster", label: "Cluster" },
  ];

  const strategy4Columns = [
    { key: "ligand_id", label: "Ligand ID" },
    { key: "smiles", label: "SMILES" },
    { key: "docking_score", label: "Docking Score" },
    {
      key: "diversity_score",
      label: "Diversity Score",
      render: (row) => {
        const value = Number(row.diversity_score);
        return Number.isFinite(value) ? value.toFixed(4) : "";
      },
    },
    {
      key: "final_score",
      label: "Final Score",
      render: (row) => {
        const value = Number(row.final_score);
        return Number.isFinite(value) ? value.toFixed(4) : "";
      },
    },
  ];

  const evaluationColumns = [
    { key: "strategy", label: "Strategy" },
    {
      key: "average_docking_score",
      label: "Average Docking Score",
      render: (row) => {
        const value = Number(row.average_docking_score);
        return Number.isFinite(value) ? value.toFixed(4) : "";
      },
    },
    {
      key: "average_similarity",
      label: "Average Similarity",
      render: (row) => {
        const value = Number(row.average_similarity);
        return Number.isFinite(value) ? value.toFixed(4) : "";
      },
    },
    { key: "redundancy_count", label: "Redundancy Count" },
  ];

  const scoreBasedRows =
    results?.strategy1_hits?.length > 0 ? results.strategy1_hits : (results?.top_hits || []);

  const dbscanCombinedRows =
    results?.strategy2_hits?.length > 0 ? results.strategy2_hits : (results?.clustered_hits || []);

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

          <h3>Docking Score Guide</h3>
          <div className="table-wrap">
            <table>
              <thead>
                <tr>
                  <th>Score</th>
                  <th>Binding strength</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td>-6</td>
                  <td>weak</td>
                </tr>
                <tr>
                  <td>-7</td>
                  <td>moderate</td>
                </tr>
                <tr>
                  <td>-8</td>
                  <td>good</td>
                </tr>
                <tr>
                  <td>-9</td>
                  <td>strong</td>
                </tr>
                <tr>
                  <td>-10 or lower</td>
                  <td>very strong</td>
                </tr>
              </tbody>
            </table>
          </div>

          <div className="download-row">
            <a href={`${API_BASE}/download-hits`} className="download-btn" target="_blank" rel="noreferrer">
              Download final_hits.csv
            </a>
            {results.strategy_evaluation_download ? (
              <a href={`${API_BASE}${results.strategy_evaluation_download}`} className="download-btn" target="_blank" rel="noreferrer">
                Download strategy_evaluation.csv
              </a>
            ) : null}
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
          {activeTab === "score_based" ? (
            <ResultsTable rows={scoreBasedRows} columns={strategyColumns} />
          ) : null}
          {activeTab === "dbscan_combined" ? (
            <ResultsTable rows={dbscanCombinedRows} columns={strategy2Columns} />
          ) : null}
          {activeTab === "final" ? <ResultsTable rows={results.final_hits || []} columns={columns} /> : null}
          {activeTab === "greedy" ? (
            <ResultsTable rows={results.strategy3_hits || []} columns={strategyColumns} />
          ) : null}
          {activeTab === "multi" ? (
            <ResultsTable rows={results.strategy4_hits || []} columns={strategy4Columns} />
          ) : null}
          {activeTab === "evaluation" ? (
            <ResultsTable rows={results.strategy_evaluation || []} columns={evaluationColumns} />
          ) : null}

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
              {results.strategy_score_plot_download ? (
                <div className="viz-card">
                  <h4>Docking Score Comparison Across Strategies</h4>
                  <img src={`${API_BASE}${results.strategy_score_plot_download}`} alt="Strategy docking score comparison" />
                </div>
              ) : null}
              {results.strategy_diversity_plot_download ? (
                <div className="viz-card">
                  <h4>Diversity Comparison Across Strategies</h4>
                  <img src={`${API_BASE}${results.strategy_diversity_plot_download}`} alt="Strategy diversity comparison" />
                </div>
              ) : null}
              {results.strategy_pca_plot_download ? (
                <div className="viz-card">
                  <h4>PCA Chemical Space for Selected Hits</h4>
                  <img src={`${API_BASE}${results.strategy_pca_plot_download}`} alt="Strategy chemical space PCA" />
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
