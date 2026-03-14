import React, { useState } from "react";

const API_BASE = "http://localhost:8000";

function UploadLigands({ onUploaded }) {
  const [file, setFile] = useState(null);
  const [status, setStatus] = useState("");

  const handleUpload = async () => {
    if (!file) {
      setStatus("Please choose a .csv file first.");
      return;
    }

    const formData = new FormData();
    formData.append("file", file);

    try {
      setStatus("Uploading ligand dataset...");
      const response = await fetch(`${API_BASE}/upload-ligands`, {
        method: "POST",
        body: formData,
      });

      if (!response.ok) {
        const error = await response.json();
        throw new Error(error.detail || "Failed to upload ligands");
      }

      const data = await response.json();
      setStatus(`Ligands uploaded. Rows detected: ${data.total_rows}`);
      onUploaded(data);
    } catch (err) {
      setStatus(`Upload failed: ${err.message}`);
    }
  };

  return (
    <section className="panel">
      <h2>Upload Ligands (.csv)</h2>
      <input type="file" accept=".csv" onChange={(e) => setFile(e.target.files?.[0] || null)} />
      <button onClick={handleUpload}>Upload Ligands</button>
      <p className="status">{status}</p>
    </section>
  );
}

export default UploadLigands;
