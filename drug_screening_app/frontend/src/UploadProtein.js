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

function UploadProtein({ onUploaded }) {
  const [file, setFile] = useState(null);
  const [status, setStatus] = useState("");
  const [uploading, setUploading] = useState(false);

  const handleUpload = async () => {
    if (uploading) {
      return;
    }

    if (!file) {
      setStatus("Please choose a .pdb file first.");
      return;
    }

    const formData = new FormData();
    formData.append("file", file);

    try {
      setUploading(true);
      setStatus("Uploading protein...");
      const response = await fetchWithTimeout(`${API_BASE}/upload-protein`, {
        method: "POST",
        body: formData,
      });

      if (!response.ok) {
        const error = await response.json();
        throw new Error(error.detail || "Failed to upload protein");
      }

      const data = await response.json();
      setStatus("Protein uploaded successfully.");
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
      <h2>Upload Protein (.pdb)</h2>
      <input type="file" accept=".pdb" onChange={(e) => setFile(e.target.files?.[0] || null)} />
      <button onClick={handleUpload} disabled={uploading}>{uploading ? "Uploading..." : "Upload Protein"}</button>
      <p className="status">{status}</p>
    </section>
  );
}

export default UploadProtein;
