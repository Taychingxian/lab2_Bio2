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
    """
    Retrieve protein structure from PDB database.
    
    Parameters:
    -----------
    protID : str
        The PDB ID of the protein (e.g., '1A2B')
    
    Returns:
    --------
    Bio.PDB.Structure.Structure
        The protein structure object
    """
    try:
        # Create a temporary directory for storing PDB files
        temp_dir = tempfile.gettempdir()
        
        # Download the PDB file
        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(protID, pdir=temp_dir, file_format='pdb')
        
        # Parse the structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(protID, pdb_file)
        
        return structure
    
    except Exception as e:
        st.error(f"Error retrieving protein structure: {e}")
        return None


def get_structure_info(prot_structure):
    """
    Calculate structural information for a protein.
    
    Parameters:
    -----------
    prot_structure : Bio.PDB.Structure.Structure
        The protein structure object
    
    Returns:
    --------
    dict
        Dictionary containing:
        - 'com': Center of mass (numpy array of [x, y, z])
        - 'Rg': Radius of gyration (float)
        - '3dview': py3Dmol view object for 3D visualization
    """
    try:
        # Extract all atoms
        atoms = [atom for atom in prot_structure.get_atoms()]
        
        # Calculate center of mass
        coords = np.array([atom.coord for atom in atoms])
        masses = np.array([atom.mass for atom in atoms])
        total_mass = np.sum(masses)
        com = np.sum(coords * masses[:, np.newaxis], axis=0) / total_mass
        
        # Calculate radius of gyration
        # Rg = sqrt(sum(mass_i * r_i^2) / total_mass)
        # where r_i is the distance from atom i to the center of mass
        distances_squared = np.sum((coords - com)**2, axis=1)
        Rg = np.sqrt(np.sum(masses * distances_squared) / total_mass)
        
        # Create 3D view using py3Dmol
        # Save structure to temporary PDB string
        temp_pdb = tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False)
        io = PDBIO()
        io.set_structure(prot_structure)
        io.save(temp_pdb.name)
        
        # Read the PDB file content
        with open(temp_pdb.name, 'r') as f:
            pdb_string = f.read()
        
        # Clean up temporary file
        os.unlink(temp_pdb.name)
        
        # Create py3Dmol view
        view = py3Dmol.view(width=800, height=600)
        view.addModel(pdb_string, 'pdb')
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        view.setBackgroundColor('white')
        view.zoomTo()
        
        # Return dictionary with all information
        return {
            'com': com,
            'Rg': Rg,
            '3dview': view
        }
    
    except Exception as e:
        st.error(f"Error calculating structure information: {e}")
        return None


def main():
    """
    Main function for the Streamlit application.
    """
    # Set page configuration
    st.set_page_config(
        page_title="Protein Structure Viewer",
        page_icon="🧬",
        layout="wide"
    )
    
    # Title and description
    st.title("🧬 Protein Structure Visualization App")
    st.markdown("""
    This application retrieves and displays structural information about proteins from the PDB database.
    Enter a valid PDB ID to view the protein's:
    - **Center of Mass (COM)**
    - **Radius of Gyration (Rg)**
    - **3D Structure Visualization**
    """)
    
    # Sidebar for input
    st.sidebar.header("Input Parameters")
    
    # Input for PDB ID
    pdb_id = st.sidebar.text_input(
        "Enter PDB ID:",
        value="",
        max_chars=4,
        help="Enter a 4-character PDB ID (e.g., 1A2B, 3NIR, 2HHB)"
    ).strip().upper()
    
    # Example PDB IDs
    st.sidebar.markdown("### Example PDB IDs:")
    st.sidebar.markdown("""
    - **1A2B** - Hemoglobin
    - **3NIR** - Insulin
    - **2HHB** - Deoxyhemoglobin
    - **1XYZ** - Example structure
    """)
    
    # Process button
    if st.sidebar.button("🔍 Analyze Protein", type="primary"):
        if len(pdb_id) == 4:
            with st.spinner(f"Retrieving protein structure for {pdb_id}..."):
                # Get protein structure
                structure = get_protein_structure(pdb_id)
                
                if structure is not None:
                    st.success(f"✅ Successfully retrieved structure for {pdb_id}")
                    
                    with st.spinner("Calculating structural properties..."):
                        # Get structure information
                        info = get_structure_info(structure)
                        
                        if info is not None:
                            # Display results in columns
                            col1, col2 = st.columns(2)
                            
                            with col1:
                                st.subheader("📊 Structural Properties")
                                
                                # Display Center of Mass
                                st.markdown("### Center of Mass (COM)")
                                st.write(f"**X:** {info['com'][0]:.3f} Å")
                                st.write(f"**Y:** {info['com'][1]:.3f} Å")
                                st.write(f"**Z:** {info['com'][2]:.3f} Å")
                                
                                st.markdown("---")
                                
                                # Display Radius of Gyration
                                st.markdown("### Radius of Gyration (Rg)")
                                st.write(f"**Rg:** {info['Rg']:.3f} Å")
                                
                                st.markdown("---")
                                
                                # Additional information
                                st.markdown("### Structure Information")
                                st.write(f"**PDB ID:** {pdb_id}")
                                
                                # Count atoms
                                atom_count = len(list(structure.get_atoms()))
                                st.write(f"**Total Atoms:** {atom_count}")
                                
                                # Count residues
                                residue_count = len(list(structure.get_residues()))
                                st.write(f"**Total Residues:** {residue_count}")
                            
                            with col2:
                                st.subheader("🧪 3D Structure Visualization")
                                # Display 3D view using HTML component
                                view_html = info['3dview']._make_html()
                                components.html(view_html, height=600, width=800, scrolling=False)
                                
                                st.info("""
                                **Interaction Tips:**
                                - 🖱️ Click and drag to rotate
                                - 📜 Scroll to zoom in/out
                                - 🔄 Right-click and drag to pan
                                """)
                else:
                    st.error(f"❌ Could not retrieve structure for {pdb_id}. Please check the PDB ID and try again.")
        else:
            st.warning("⚠️ Please enter a valid 4-character PDB ID.")
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center'>
        <p>Developed by Tay Ching Xian | Bio2 Lab 2</p>
        <p>Data source: <a href='https://www.rcsb.org/' target='_blank'>RCSB Protein Data Bank</a></p>
    </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
