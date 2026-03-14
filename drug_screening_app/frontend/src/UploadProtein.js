import React, { useState } from "react";

const API_BASE = "http://localhost:8000";

function UploadProtein({ onUploaded }) {
  const [file, setFile] = useState(null);
  const [status, setStatus] = useState("");

  const handleUpload = async () => {
    if (!file) {
      setStatus("Please choose a .pdb file first.");
      return;
    }

    const formData = new FormData();
    formData.append("file", file);

    try {
      setStatus("Uploading protein...");
      const response = await fetch(`${API_BASE}/upload-protein`, {
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
      setStatus(`Upload failed: ${err.message}`);
    }
  };

  return (
    <section className="panel">
      <h2>Upload Protein (.pdb)</h2>
      <input type="file" accept=".pdb" onChange={(e) => setFile(e.target.files?.[0] || null)} />
      <button onClick={handleUpload}>Upload Protein</button>
      <p className="status">{status}</p>
    </section>
  );
}

export default UploadProtein;
