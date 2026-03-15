import React, { useState } from "react";

const API_BASE = "http://127.0.0.1:8000";
const REQUEST_TIMEOUT_MS = 30000;

async function fetchWithTimeout(url, options = {}, timeoutMs = REQUEST_TIMEOUT_MS) {
  const controller = new AbortController();
  const timer = setTimeout(() => controller.abort(), timeoutMs);

  try {
    return await fetch(url, { ...options, signal: controller.signal });
  } finally {
    clearTimeout(timer);
  }
}

function UploadLigands({ onUploaded }) {
  const [file, setFile] = useState(null);
  const [status, setStatus] = useState("");
  const [uploading, setUploading] = useState(false);

  const handleUpload = async () => {
    if (uploading) {
      return;
    }

    if (!file) {
      setStatus("Please choose a .csv file first.");
      return;
    }

    const formData = new FormData();
    formData.append("file", file);

    try {
      setUploading(true);
      setStatus("Uploading ligand dataset...");
      const response = await fetchWithTimeout(`${API_BASE}/upload-ligands`, {
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
      const isTimeout = err.name === "AbortError";
      setStatus(isTimeout ? "Upload timed out. Check backend and try again." : `Upload failed: ${err.message}`);
    } finally {
      setUploading(false);
    }
  };

  return (
    <section className="panel">
      <h2>Upload Ligands (.csv)</h2>
      <input type="file" accept=".csv" onChange={(e) => setFile(e.target.files?.[0] || null)} />
      <button onClick={handleUpload} disabled={uploading}>{uploading ? "Uploading..." : "Upload Ligands"}</button>
      <p className="status">{status}</p>
    </section>
  );
}

export default UploadLigands;
