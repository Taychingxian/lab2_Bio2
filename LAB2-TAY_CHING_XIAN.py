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
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    # Custom CSS for better UI
    st.markdown("""
        <style>
        /* Main container styling */
        .main {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 0;
        }
        
        /* Card styling */
        .stMarkdown {
            background-color: rgba(255, 255, 255, 0.95);
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }
        
        /* Header styling */
        h1 {
            color: #1e3a8a;
            font-weight: 700;
            text-align: center;
            padding: 20px 0;
            text-shadow: 2px 2px 4px rgba(0, 0, 0, 0.1);
        }
        
        h2, h3 {
            color: #2563eb;
            font-weight: 600;
        }
        
        /* Sidebar styling */
        [data-testid="stSidebar"] {
            background: linear-gradient(180deg, #1e3a8a 0%, #3b82f6 100%);
        }
        
        [data-testid="stSidebar"] [data-testid="stMarkdownContainer"] {
            color: white;
        }
        
        [data-testid="stSidebar"] h1, 
        [data-testid="stSidebar"] h2, 
        [data-testid="stSidebar"] h3 {
            color: white !important;
        }
        
        [data-testid="stSidebar"] label {
            color: white !important;
            font-weight: 500;
        }
        
        /* Input field styling */
        .stTextInput input {
            border: 2px solid #3b82f6;
            border-radius: 8px;
            padding: 10px;
            font-size: 16px;
            transition: all 0.3s ease;
        }
        
        .stTextInput input:focus {
            border-color: #2563eb;
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.1);
        }
        
        /* Button styling */
        .stButton button {
            background: linear-gradient(135deg, #10b981 0%, #059669 100%);
            color: white;
            border: none;
            border-radius: 8px;
            padding: 12px 24px;
            font-size: 16px;
            font-weight: 600;
            transition: all 0.3s ease;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
            width: 100%;
        }
        
        .stButton button:hover {
            transform: translateY(-2px);
            box-shadow: 0 6px 12px rgba(0, 0, 0, 0.15);
        }
        
        /* Success/Error message styling */
        .stSuccess {
            background-color: #d1fae5;
            border-left: 4px solid #10b981;
            border-radius: 8px;
            padding: 15px;
            animation: slideIn 0.5s ease;
        }
        
        .stError {
            background-color: #fee2e2;
            border-left: 4px solid #ef4444;
            border-radius: 8px;
            padding: 15px;
            animation: slideIn 0.5s ease;
        }
        
        .stWarning {
            background-color: #fef3c7;
            border-left: 4px solid #f59e0b;
            border-radius: 8px;
            padding: 15px;
            animation: slideIn 0.5s ease;
        }
        
        .stInfo {
            background-color: #dbeafe;
            border-left: 4px solid #3b82f6;
            border-radius: 8px;
            padding: 15px;
        }
        
        /* Metric styling */
        [data-testid="stMetricValue"] {
            font-size: 28px;
            font-weight: 700;
            color: #1e3a8a;
        }
        
        [data-testid="stMetricLabel"] {
            font-size: 14px;
            font-weight: 500;
            color: #64748b;
        }
        
        /* Animations */
        @keyframes slideIn {
            from {
                opacity: 0;
                transform: translateX(-20px);
            }
            to {
                opacity: 1;
                transform: translateX(0);
            }
        }
        
        @keyframes fadeIn {
            from { opacity: 0; }
            to { opacity: 1; }
        }
        
        /* Property cards */
        .property-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
            margin: 10px 0;
            animation: fadeIn 0.5s ease;
        }
        
        .property-value {
            font-size: 24px;
            font-weight: 700;
            margin: 5px 0;
        }
        
        .property-label {
            font-size: 14px;
            opacity: 0.9;
            text-transform: uppercase;
            letter-spacing: 1px;
        }
        
        /* 3D viewer container */
        .viewer-container {
            background: white;
            border-radius: 12px;
            padding: 15px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
            animation: fadeIn 0.7s ease;
        }
        
        /* Footer styling */
        .footer {
            text-align: center;
            padding: 20px;
            color: #64748b;
            font-size: 14px;
            margin-top: 40px;
        }
        
        /* Divider */
        hr {
            border: none;
            height: 2px;
            background: linear-gradient(90deg, transparent, #3b82f6, transparent);
            margin: 30px 0;
        }
        </style>
    """, unsafe_allow_html=True)
    
    # Hero Section with gradient background
    st.markdown("""
        <div style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                    padding: 40px; border-radius: 15px; margin-bottom: 30px; 
                    box-shadow: 0 8px 16px rgba(0, 0, 0, 0.2);'>
            <h1 style='color: white; font-size: 48px; margin: 0; text-shadow: 2px 2px 4px rgba(0,0,0,0.3);'>
                🧬 Protein Structure Viewer
            </h1>
            <p style='color: rgba(255,255,255,0.95); font-size: 18px; text-align: center; margin-top: 15px;'>
                Explore the molecular architecture of life
            </p>
        </div>
    """, unsafe_allow_html=True)
    
    # Introduction with icons
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("""
            <div style='text-align: center; padding: 20px; background: linear-gradient(135deg, #3b82f6 0%, #2563eb 100%); 
                        border-radius: 12px; color: white; box-shadow: 0 4px 8px rgba(0,0,0,0.1);'>
                <h2 style='color: white; font-size: 40px; margin: 10px 0;'>📍</h2>
                <h3 style='color: white; font-size: 18px; margin: 10px 0;'>Center of Mass</h3>
                <p style='font-size: 14px; opacity: 0.9;'>Geometric center weighted by atomic masses</p>
            </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
            <div style='text-align: center; padding: 20px; background: linear-gradient(135deg, #10b981 0%, #059669 100%); 
                        border-radius: 12px; color: white; box-shadow: 0 4px 8px rgba(0,0,0,0.1);'>
                <h2 style='color: white; font-size: 40px; margin: 10px 0;'>📏</h2>
                <h3 style='color: white; font-size: 18px; margin: 10px 0;'>Radius of Gyration</h3>
                <p style='font-size: 14px; opacity: 0.9;'>Measure of protein compactness</p>
            </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown("""
            <div style='text-align: center; padding: 20px; background: linear-gradient(135deg, #f59e0b 0%, #d97706 100%); 
                        border-radius: 12px; color: white; box-shadow: 0 4px 8px rgba(0,0,0,0.1);'>
                <h2 style='color: white; font-size: 40px; margin: 10px 0;'>🔬</h2>
                <h3 style='color: white; font-size: 18px; margin: 10px 0;'>3D Visualization</h3>
                <p style='font-size: 14px; opacity: 0.9;'>Interactive molecular structure</p>
            </div>
        """, unsafe_allow_html=True)
    
    st.markdown("<br>", unsafe_allow_html=True)
    
    # Sidebar for input
    st.sidebar.markdown("""
        <div style='text-align: center; padding: 20px 0;'>
            <h2 style='color: white; font-size: 24px; margin: 0;'>⚙️ Control Panel</h2>
        </div>
    """, unsafe_allow_html=True)
    
    # Input for PDB ID
    st.sidebar.markdown("<br>", unsafe_allow_html=True)
    pdb_id = st.sidebar.text_input(
        "🔍 Enter PDB ID:",
        value="",
        max_chars=4,
        help="Enter a 4-character PDB ID (e.g., 1A2B, 3NIR, 2HHB)",
        placeholder="e.g., 1A2B"
    ).strip().upper()
    
    # Example PDB IDs with enhanced styling
    st.sidebar.markdown("<br>", unsafe_allow_html=True)
    st.sidebar.markdown("""
        <div style='background: rgba(255,255,255,0.1); padding: 15px; border-radius: 10px; margin: 10px 0;'>
            <h3 style='color: white; font-size: 18px; margin-bottom: 10px;'>💡 Example Proteins</h3>
    """, unsafe_allow_html=True)
    
    examples = {
        "1A2B": "Hemoglobin",
        "3NIR": "Insulin",
        "2HHB": "Deoxyhemoglobin",
        "1MBN": "Myoglobin"
    }
    
    for pdb, name in examples.items():
        st.sidebar.markdown(f"""
            <div style='padding: 8px; margin: 5px 0; background: rgba(255,255,255,0.15); 
                        border-radius: 6px; border-left: 3px solid #10b981;'>
                <strong style='color: #fff;'>{pdb}</strong> 
                <span style='color: rgba(255,255,255,0.8);'>- {name}</span>
            </div>
        """, unsafe_allow_html=True)
    
    st.sidebar.markdown("</div>", unsafe_allow_html=True)
    
    # Process button with icon
    st.sidebar.markdown("<br>", unsafe_allow_html=True)
    analyze_button = st.sidebar.button("🔬 Analyze Protein Structure", type="primary", use_container_width=True)
    
    # Additional info in sidebar
    st.sidebar.markdown("<br><br>", unsafe_allow_html=True)
    st.sidebar.markdown("""
        <div style='background: rgba(255,255,255,0.1); padding: 15px; border-radius: 10px;'>
            <p style='color: white; font-size: 12px; margin: 0; line-height: 1.6;'>
                <strong>ℹ️ About:</strong><br>
                Data from RCSB Protein Data Bank<br>
                Developed by Tay Ching Xian<br>
                Bio2 Lab 2 - 2025
            </p>
        </div>
    """, unsafe_allow_html=True)
    
    if analyze_button:
        if len(pdb_id) == 4:
            with st.spinner(f"🔄 Retrieving protein structure for **{pdb_id}**..."):
                # Get protein structure
                structure = get_protein_structure(pdb_id)
                
                if structure is not None:
                    st.markdown(f"""
                        <div style='background: linear-gradient(135deg, #10b981 0%, #059669 100%); 
                                    padding: 20px; border-radius: 12px; margin: 20px 0; 
                                    box-shadow: 0 4px 8px rgba(0,0,0,0.1); animation: slideIn 0.5s ease;'>
                            <h3 style='color: white; margin: 0; font-size: 20px;'>
                                ✅ Successfully Retrieved: {pdb_id}
                            </h3>
                            <p style='color: rgba(255,255,255,0.9); margin: 10px 0 0 0;'>
                                Processing structural data...
                            </p>
                        </div>
                    """, unsafe_allow_html=True)
                    
                    with st.spinner("⚡ Calculating structural properties..."):
                        # Get structure information
                        info = get_structure_info(structure)
                        
                        if info is not None:
                            # Display results
                            st.markdown("<br>", unsafe_allow_html=True)
                            
                            # Main content area with two columns
                            col1, col2 = st.columns([1, 1.5], gap="large")
                            
                            with col1:
                                st.markdown("""
                                    <div style='background: white; padding: 25px; border-radius: 15px; 
                                                box-shadow: 0 4px 12px rgba(0,0,0,0.1); animation: fadeIn 0.5s ease;'>
                                        <h2 style='color: #1e3a8a; margin-top: 0; display: flex; align-items: center;'>
                                            📊 Structural Properties
                                        </h2>
                                    </div>
                                """, unsafe_allow_html=True)
                                
                                st.markdown("<br>", unsafe_allow_html=True)
                                
                                # Center of Mass Card
                                st.markdown("""
                                    <div class='property-card'>
                                        <div class='property-label'>📍 Center of Mass (COM)</div>
                                """, unsafe_allow_html=True)
                                
                                # Display COM values in metrics
                                col_x, col_y, col_z = st.columns(3)
                                with col_x:
                                    st.metric("X Coordinate", f"{info['com'][0]:.2f} Å")
                                with col_y:
                                    st.metric("Y Coordinate", f"{info['com'][1]:.2f} Å")
                                with col_z:
                                    st.metric("Z Coordinate", f"{info['com'][2]:.2f} Å")
                                
                                st.markdown("</div>", unsafe_allow_html=True)
                                
                                st.markdown("<br>", unsafe_allow_html=True)
                                
                                # Radius of Gyration Card
                                st.markdown(f"""
                                    <div class='property-card'>
                                        <div class='property-label'>📏 Radius of Gyration</div>
                                        <div class='property-value'>{info['Rg']:.2f} Å</div>
                                        <p style='font-size: 13px; margin-top: 10px; opacity: 0.9;'>
                                            Measure of structural compactness
                                        </p>
                                    </div>
                                """, unsafe_allow_html=True)
                                
                                st.markdown("<br>", unsafe_allow_html=True)
                                
                                # Additional Information Card
                                atom_count = len(list(structure.get_atoms()))
                                residue_count = len(list(structure.get_residues()))
                                
                                st.markdown("""
                                    <div style='background: linear-gradient(135deg, #f59e0b 0%, #d97706 100%); 
                                                color: white; padding: 20px; border-radius: 12px; 
                                                box-shadow: 0 4px 12px rgba(0,0,0,0.15); animation: fadeIn 0.6s ease;'>
                                        <div class='property-label'>📋 Structure Details</div>
                                """, unsafe_allow_html=True)
                                
                                info_col1, info_col2 = st.columns(2)
                                with info_col1:
                                    st.markdown(f"""
                                        <div style='text-align: center; padding: 10px;'>
                                            <div style='font-size: 28px; font-weight: 700;'>{atom_count:,}</div>
                                            <div style='font-size: 12px; opacity: 0.9;'>ATOMS</div>
                                        </div>
                                    """, unsafe_allow_html=True)
                                with info_col2:
                                    st.markdown(f"""
                                        <div style='text-align: center; padding: 10px;'>
                                            <div style='font-size: 28px; font-weight: 700;'>{residue_count:,}</div>
                                            <div style='font-size: 12px; opacity: 0.9;'>RESIDUES</div>
                                        </div>
                                    """, unsafe_allow_html=True)
                                
                                st.markdown(f"""
                                        <div style='text-align: center; margin-top: 15px; padding-top: 15px; 
                                                    border-top: 1px solid rgba(255,255,255,0.3);'>
                                            <strong>PDB ID:</strong> {pdb_id}
                                        </div>
                                    </div>
                                """, unsafe_allow_html=True)
                            
                            with col2:
                                st.markdown("""
                                    <div class='viewer-container'>
                                        <h2 style='color: #1e3a8a; margin-top: 0; display: flex; align-items: center;'>
                                            🔬 3D Structure Visualization
                                        </h2>
                                """, unsafe_allow_html=True)
                                
                                # Display 3D view using HTML component
                                view_html = info['3dview']._make_html()
                                components.html(view_html, height=550, scrolling=False)
                                
                                st.markdown("""
                                    <div style='background: linear-gradient(135deg, #3b82f6 0%, #2563eb 100%); 
                                                color: white; padding: 15px; border-radius: 10px; margin-top: 15px;'>
                                        <strong>💡 Interaction Tips:</strong>
                                        <ul style='margin: 10px 0 0 0; padding-left: 20px; line-height: 1.8;'>
                                            <li>🖱️ <strong>Rotate:</strong> Click and drag</li>
                                            <li>🔍 <strong>Zoom:</strong> Scroll wheel</li>
                                            <li>↔️ <strong>Pan:</strong> Right-click and drag</li>
                                        </ul>
                                    </div>
                                    </div>
                                """, unsafe_allow_html=True)
                else:
                    st.markdown(f"""
                        <div style='background: linear-gradient(135deg, #ef4444 0%, #dc2626 100%); 
                                    padding: 20px; border-radius: 12px; margin: 20px 0; 
                                    box-shadow: 0 4px 8px rgba(0,0,0,0.1); animation: slideIn 0.5s ease;'>
                            <h3 style='color: white; margin: 0; font-size: 20px;'>
                                ❌ Error: Could Not Retrieve Structure
                            </h3>
                            <p style='color: rgba(255,255,255,0.9); margin: 10px 0 0 0;'>
                                Please verify that <strong>{pdb_id}</strong> is a valid PDB ID and try again.
                            </p>
                        </div>
                    """, unsafe_allow_html=True)
        else:
            st.markdown("""
                <div style='background: linear-gradient(135deg, #f59e0b 0%, #d97706 100%); 
                            padding: 20px; border-radius: 12px; margin: 20px 0; 
                            box-shadow: 0 4px 8px rgba(0,0,0,0.1); animation: slideIn 0.5s ease;'>
                    <h3 style='color: white; margin: 0; font-size: 20px;'>
                        ⚠️ Invalid Input
                    </h3>
                    <p style='color: rgba(255,255,255,0.9); margin: 10px 0 0 0;'>
                        Please enter a valid 4-character PDB ID (e.g., 1A2B, 3NIR, 2HHB)
                    </p>
                </div>
            """, unsafe_allow_html=True)
    
    # Footer
    st.markdown("<br><br>", unsafe_allow_html=True)
    st.markdown("""
        <div style='background: white; padding: 30px; border-radius: 15px; 
                    box-shadow: 0 4px 12px rgba(0,0,0,0.1); text-align: center; margin-top: 40px;'>
            <div style='display: flex; justify-content: center; align-items: center; gap: 30px; flex-wrap: wrap;'>
                <div>
                    <div style='color: #64748b; font-size: 14px; margin-bottom: 5px;'>Developed by</div>
                    <div style='color: #1e3a8a; font-weight: 600; font-size: 16px;'>Tay Ching Xian</div>
                </div>
                <div style='width: 2px; height: 40px; background: #e2e8f0;'></div>
                <div>
                    <div style='color: #64748b; font-size: 14px; margin-bottom: 5px;'>Course</div>
                    <div style='color: #1e3a8a; font-weight: 600; font-size: 16px;'>Bio2 Lab 2</div>
                </div>
                <div style='width: 2px; height: 40px; background: #e2e8f0;'></div>
                <div>
                    <div style='color: #64748b; font-size: 14px; margin-bottom: 5px;'>Data Source</div>
                    <div style='color: #1e3a8a; font-weight: 600; font-size: 16px;'>
                        <a href='https://www.rcsb.org/' target='_blank' 
                           style='color: #3b82f6; text-decoration: none;'>RCSB PDB</a>
                    </div>
                </div>
            </div>
            <div style='margin-top: 20px; padding-top: 20px; border-top: 2px solid #e2e8f0;'>
                <p style='color: #94a3b8; font-size: 13px; margin: 0;'>
                    © 2025 | Protein Structure Visualization Tool | Educational Use Only
                </p>
            </div>
        </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
