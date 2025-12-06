# lab2_Bio2
# 🧬 Protein Structure Viewer

## 📖 Description
This **Protein Structure Viewer** is an interactive web application built with [Streamlit](https://streamlit.io/). It allows users to retrieve biological macromolecular structures directly from the **RCSB Protein Data Bank (PDB)**, visualize them in 3D, and calculate key structural properties in real-time.

The tool is designed for students, researchers, and bioinformaticians who need a quick, lightweight interface to assess the geometry and compactness of protein structures.

## ✨ Key Features

* **⚡ Real-time Data Retrieval:** Fetches structure files (`.pdb`) instantly using a 4-character PDB ID (e.g., `1A2B`, `4HHB`).
* **🔬 Interactive 3D Visualization:**
    * Rendered using **py3Dmol**.
    * Supports Zooming, Panning, and Rotating.
    * Visualizes protein backbone in "Cartoon" style with Spectrum coloring.
* **📊 Structural Metrics:** Automatically calculates:
    * **Center of Mass (COM):** The spatial average of the protein's mass.
    * **Radius of Gyration ($R_g$):** A measure of the protein's compactness.
    * **Composition:** Total count of Atoms and Residues.

## 🛠️ Tech Stack
* **Python 3.x**
* **Streamlit:** For the web interface.
* **BioPython (`Bio.PDB`):** For parsing PDB files and handling atomic data.
* **py3Dmol:** For embedding the 3D molecular viewer.
* **NumPy:** For vector mathematics and matrix calculations.
