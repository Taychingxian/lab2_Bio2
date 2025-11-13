"""
Protein Structure Visualization App
Author: Tay Ching Xian
Description: A Streamlit application that displays protein structural information
             including center of mass, radius of gyration, and 3D visualization.
"""

import streamlit as st
import py3Dmol
import streamlit.components.v1 as components
from Bio.PDB import PDBParser, Structure, PDBIO
from Bio.PDB.PDBList import PDBList
import numpy as np
import os
import tempfile


def get_protein_structure(protID):
    try:
        temp_dir = tempfile.gettempdir()
        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(protID, pdir=temp_dir, file_format='pdb')
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(protID, pdb_file)
        return structure
    except Exception as e:
        st.error(f"Error retrieving protein structure: {e}")
        return None


def get_structure_info(prot_structure):
    try:
        atoms = [atom for atom in prot_structure.get_atoms()]
        coords = np.array([atom.coord for atom in atoms])
        masses = np.array([atom.mass for atom in atoms])
        total_mass = np.sum(masses)
        com = np.sum(coords * masses[:, np.newaxis], axis=0) / total_mass
        distances_squared = np.sum((coords - com)**2, axis=1)
        Rg = np.sqrt(np.sum(masses * distances_squared) / total_mass)
        
        temp_pdb = tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False)
        io = PDBIO()
        io.set_structure(prot_structure)
        io.save(temp_pdb.name)
        
        with open(temp_pdb.name, 'r') as f:
            pdb_string = f.read()
        
        os.unlink(temp_pdb.name)
        
        view = py3Dmol.view(width=800, height=600)
        view.addModel(pdb_string, 'pdb')
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        view.setBackgroundColor('#1e293b')
        view.zoomTo()
        view._model_data = pdb_string
        
        return {
            'com': com,
            'Rg': Rg,
            '3dview': view
        }
    except Exception as e:
        st.error(f"Error calculating structure information: {e}")
        return None


def main():
    st.set_page_config(
        page_title="Protein Structure Viewer",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    st.markdown("""
        <style>
        [data-testid="stSidebar"] {
            background-color: #1e293b;
        }
        [data-testid="stSidebar"] label {
            color: white !important;
            font-weight: 500;
        }
        .stButton button {
            background-color: #0ea5e9;
            color: white;
            border: none;
            border-radius: 8px;
            padding: 12px 24px;
            font-size: 16px;
            font-weight: 600;
            width: 100%;
        }
        .stButton button:hover {
            background-color: #0284c7;
        }
        [data-testid="stMetricValue"] {
            font-size: 24px;
            color: #0ea5e9;
        }
        hr {
            margin: 20px 0;
        }
        </style>
    """, unsafe_allow_html=True)
    
    st.title("🧬 Protein Structure Viewer")
    st.markdown("Analyze protein structures from the PDB database")
    st.markdown("---")
    
    st.sidebar.header("Input Parameters")
    
    pdb_id = st.sidebar.text_input(
        "Enter PDB ID:",
        value="",
        max_chars=4,
        help="Enter a 4-character PDB ID (e.g., 1A2B, 3NIR, 2HHB)",
        placeholder="e.g., 1A2B"
    ).strip().upper()
    
    st.sidebar.markdown("### Examples:")
    st.sidebar.markdown("- **1A2B** - Hemoglobin")
    st.sidebar.markdown("- **3NIR** - Insulin")
    st.sidebar.markdown("- **2HHB** - Deoxyhemoglobin")
    st.sidebar.markdown("- **1MBN** - Myoglobin")
    
    analyze_button = st.sidebar.button("🔬 Analyze Protein", type="primary", use_container_width=True)
    
    st.sidebar.markdown("---")
    st.sidebar.caption("Developed by Tay Ching Xian")
    
    if analyze_button:
        if len(pdb_id) == 4:
            with st.spinner(f"Retrieving protein structure for {pdb_id}..."):
                structure = get_protein_structure(pdb_id)
                
                if structure is not None:
                    st.success(f"✅ Successfully retrieved structure for {pdb_id}")
                    
                    with st.spinner("Calculating structural properties..."):
                        info = get_structure_info(structure)
                        
                        if info is not None:
                            col1, col2 = st.columns([1, 1.5], gap="large")
                            
                            with col1:
                                st.subheader("📊 Structural Properties")
                                
                                st.markdown("#### Center of Mass (COM)")
                                col_x, col_y, col_z = st.columns(3)
                                with col_x:
                                    st.metric("X", f"{info['com'][0]:.2f} Å")
                                with col_y:
                                    st.metric("Y", f"{info['com'][1]:.2f} Å")
                                with col_z:
                                    st.metric("Z", f"{info['com'][2]:.2f} Å")
                                
                                st.markdown("---")
                                
                                st.markdown("#### Radius of Gyration")
                                st.metric("Rg", f"{info['Rg']:.2f} Å")
                                
                                st.markdown("---")
                                
                                st.markdown("#### Structure Details")
                                atom_count = len(list(structure.get_atoms()))
                                residue_count = len(list(structure.get_residues()))
                                
                                info_col1, info_col2 = st.columns(2)
                                with info_col1:
                                    st.metric("Atoms", f"{atom_count:,}")
                                with info_col2:
                                    st.metric("Residues", f"{residue_count:,}")
                                
                                st.info(f"**PDB ID:** {pdb_id}")
                            
                            with col2:
                                st.subheader("🔬 3D Structure Visualization")
                                
                                view_html = info['3dview']._make_html()
                                
                                st.markdown("""
                                    <div style='border: 2px solid #e2e8f0; border-radius: 10px; 
                                                padding: 10px; background-color: #f8fafc;'>
                                """, unsafe_allow_html=True)
                                
                                components.html(view_html, height=520, scrolling=False)
                                
                                st.markdown("</div>", unsafe_allow_html=True)
                                
                                st.markdown("<br>", unsafe_allow_html=True)
                                
                                control_col1, control_col2 = st.columns(2)
                                
                                with control_col1:
                                    st.markdown("""
                                        <div style='background: linear-gradient(135deg, #0ea5e9 0%, #0284c7 100%); 
                                                    color: white; padding: 15px; border-radius: 8px;'>
                                            <strong>🖱️ Mouse Controls</strong>
                                            <ul style='margin: 10px 0 0 0; padding-left: 20px; font-size: 14px;'>
                                                <li>Left Click + Drag: Rotate</li>
                                                <li>Right Click + Drag: Pan</li>
                                                <li>Scroll Wheel: Zoom</li>
                                            </ul>
                                        </div>
                                    """, unsafe_allow_html=True)
                                
                                with control_col2:
                                    st.markdown("""
                                        <div style='background: linear-gradient(135deg, #8b5cf6 0%, #7c3aed 100%); 
                                                    color: white; padding: 15px; border-radius: 8px;'>
                                            <strong>ℹ️ View Options</strong>
                                            <ul style='margin: 10px 0 0 0; padding-left: 20px; font-size: 14px;'>
                                                <li>Cartoon: Secondary structure</li>
                                                <li>Stick: Atomic bonds</li>
                                                <li>Surface: Molecular surface</li>
                                            </ul>
                                        </div>
                                    """, unsafe_allow_html=True)
                else:
                    st.error(f"❌ Could not retrieve structure for {pdb_id}. Please check the PDB ID.")
        else:
            st.warning("⚠️ Please enter a valid 4-character PDB ID.")
    
    st.markdown("---")
    st.markdown("""
        <div style='text-align: center; color: #64748b;'>
            <p>Data from <a href='https://www.rcsb.org/' target='_blank'>RCSB Protein Data Bank</a></p>
        </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
