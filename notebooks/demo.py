#!/usr/bin/env python3
"""
PEMFC Bio-Inspired CFD Study - Demo Script
==========================================

This script demonstrates the key features of the PEMFC CFD study:
1. Geometry generation and meshing
2. Physics setup and material properties
3. Solver configuration
4. Results analysis and visualization

Prerequisites:
Make sure you have installed all required packages:
pip install -r requirements.txt
"""

# Import required packages
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Add src to path
sys.path.append(os.path.join(os.path.dirname('.'), '..'))

# Import PEMFC modules
from geometry.flowfield import PEMFCGeometry
from src.physics import PEMFCPhysics, MaterialProperties, get_default_properties
from src.solver import PEMFCSolver

# Set plotting style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

print("PEMFC CFD Study Demo - All packages imported successfully!")

def main():
    """Main demo function."""
    
    print("\n" + "="*60)
    print("PEMFC BIO-INSPIRED CFD STUDY DEMO")
    print("="*60)
    
    # 1. Geometry Generation
    print("\n1. CREATING GEOMETRY")
    print("-" * 30)
    
    pemfc_geom = PEMFCGeometry(
        channel_depth=1.0e-3,           # 1.0 mm
        channel_width_inlet=1.0e-3,     # 1.0 mm
        channel_width_outlet=0.85e-3,   # 0.85 mm (15% taper)
        rib_width=1.0e-3,               # 1.0 mm
        gdl_thickness=2.5e-4,           # 0.25 mm
        membrane_thickness=7.5e-5,      # 0.075 mm
        num_channels=8,                  # 8 interdigitated channels
        channel_length=20.0e-3,         # 20 mm
        manifold_levels=2                # Two-level branching
    )
    
    print("Geometry parameters:")
    print(f"  Total width: {pemfc_geom.total_width*1000:.1f} mm")
    print(f"  Channel taper: {(1 - pemfc_geom.channel_width_outlet/pemfc_geom.channel_width_inlet)*100:.1f}%")
    print(f"  Active area: {pemfc_geom.total_width * pemfc_geom.channel_length * 1e4:.1f} cm²")
    
    # 2. Mesh Generation
    print("\n2. GENERATING MESHES")
    print("-" * 30)
    
    mesh_files = {}
    
    for mesh_size in ["coarse", "medium", "fine"]:
        print(f"\nGenerating {mesh_size} mesh...")
        
        try:
            # Generate mesh
            mesh_path = pemfc_geom.generate_mesh(mesh_size)
            fenics_path = pemfc_geom.export_to_fenics(mesh_path)
            
            mesh_files[mesh_size] = {
                "gmsh": mesh_path,
                "fenics": fenics_path
            }
            
            print(f"  ✓ {mesh_size} mesh generated successfully")
            print(f"  Gmsh file: {mesh_path}")
            print(f"  FEniCS file: {fenics_path}")
            
        except Exception as e:
            print(f"  ✗ Failed to generate {mesh_size} mesh: {e}")
            continue
    
    print(f"\nMesh generation complete! Generated {len(mesh_files)} meshes.")
    
    # 3. Physics Setup
    print("\n3. SETTING UP PHYSICS")
    print("-" * 30)
    
    props = get_default_properties()
    physics = PEMFCPhysics(props)
    
    print("Physics module configured with:")
    print(f"  Air density: {props.rho_air:.3f} kg/m³")
    print(f"  Air viscosity: {props.mu_air:.2e} Pa·s")
    print(f"  GDL porosity: {props.gdl_porosity:.1f}")
    print(f"  GDL permeability: {props.gdl_permeability:.1e} m²")
    print(f"  Exchange current density: {props.exchange_current_density:.1e} A/m²")
    
    # 4. Operating Conditions Analysis
    print("\n4. ANALYZING OPERATING CONDITIONS")
    print("-" * 30)
    
    current_densities = [0.3, 0.8, 1.2]  # A/cm²
    stoichiometries = [1.5, 2.5]
    active_area = 25.0e-4  # m²
    
    # Constants
    F = 96485.0  # C/mol
    R = 8.314    # J/(mol·K)
    T = 353.0    # K
    P = 1.013e5  # Pa
    x_o2_in = 0.21  # O₂ mole fraction
    
    # Calculate operating conditions
    operating_conditions = []
    
    for j in current_densities:
        for lambda_stoich in stoichiometries:
            # Convert current density to A/m²
            j_am2 = j * 1e4
            
            # Calculate O₂ consumption rate
            n_o2_cons = j_am2 * active_area / (4 * F)
            
            # Calculate inlet O₂ flow rate
            n_o2_in = lambda_stoich * n_o2_cons
            
            # Calculate total inlet flow rate
            n_tot_in = n_o2_in / x_o2_in
            
            # Calculate volumetric flow rate
            Q = n_tot_in * R * T / P
            
            # Calculate bulk velocity (assuming inlet area)
            inlet_area = 25.0e-6  # m² (approximate)
            U_bulk = Q / inlet_area
            
            # Calculate Reynolds number
            Re = props.rho_air * U_bulk * 0.001 / props.mu_air
            
            operating_conditions.append({
                'current_density': j,
                'stoichiometry': lambda_stoich,
                'o2_consumption_rate': n_o2_cons,
                'o2_inlet_rate': n_o2_in,
                'total_inlet_rate': n_tot_in,
                'volumetric_flow': Q,
                'bulk_velocity': U_bulk,
                'reynolds_number': Re
            })
    
    # Create DataFrame for analysis
    op_df = pd.DataFrame(operating_conditions)
    print("Operating Conditions Matrix:")
    print(op_df.round(6))
    
    # 5. CFD Simulation Setup
    print("\n5. SETTING UP CFD SOLVER")
    print("-" * 30)
    
    if "medium" in mesh_files:
        mesh_file = mesh_files["medium"]["fenics"]
        mesh_refinement = "medium"
        
        print(f"Setting up solver for {mesh_refinement} mesh...")
        print(f"Mesh file: {mesh_file}")
        
        try:
            # Initialize solver
            solver = PEMFCSolver(mesh_file, physics, mesh_refinement)
            print("✓ Solver initialized successfully")
            print(f"  Mesh cells: {solver.mesh.topology.index_map(solver.mesh.topology.dim).size_global:,}")
            print(f"  Mesh vertices: {solver.mesh.topology.index_map(0).size_global:,}")
            
        except Exception as e:
            print(f"✗ Failed to initialize solver: {e}")
            solver = None
    else:
        print("No medium mesh available for solver setup")
        solver = None
    
    # 6. Results Analysis
    print("\n6. ANALYZING RESULTS")
    print("-" * 30)
    
    print("Operating Conditions Summary:")
    print("=" * 50)
    
    # Flow regime analysis
    print("\nFlow Regime Analysis:")
    for _, row in op_df.iterrows():
        regime = "Laminar" if row['reynolds_number'] < 2300 else "Transitional"
        print(f"  j={row['current_density']} A/cm², λ={row['stoichiometry']}: "
              f"Re={row['reynolds_number']:.1f} ({regime})")
    
    # Mass flow analysis
    print("\nMass Flow Analysis:")
    for _, row in op_df.iterrows():
        print(f"  j={row['current_density']} A/cm², λ={row['stoichiometry']}: "
              f"ṁ_O2={row['o2_inlet_rate']*32e-3:.3f} g/s")
    
    # Power analysis
    print("\nPower Analysis:")
    for _, row in op_df.iterrows():
        power = row['current_density'] * active_area * 1.23  # Assuming 1.23V
        print(f"  j={row['current_density']} A/cm²: P={power:.2f} W")
    
    # 7. Visualization and Post-Processing
    print("\n7. CREATING VISUALIZATIONS")
    print("-" * 30)
    
    # Create comprehensive visualization
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('PEMFC Operating Conditions Analysis', fontsize=16, fontweight='bold')
    
    # 1. Reynolds number vs current density
    ax1 = axes[0, 0]
    for lambda_stoich in stoichiometries:
        subset = op_df[op_df['stoichiometry'] == lambda_stoich]
        ax1.plot(subset['current_density'], subset['reynolds_number'], 
                  'o-', linewidth=2, markersize=8, label=f'λ = {lambda_stoich}')
    ax1.set_xlabel('Current Density [A/cm²]')
    ax1.set_ylabel('Reynolds Number')
    ax1.set_title('Flow Regime')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # 2. Volumetric flow rate vs current density
    ax2 = axes[0, 1]
    for lambda_stoich in stoichiometries:
        subset = op_df[op_df['stoichiometry'] == lambda_stoich]
        ax2.plot(subset['current_density'], subset['volumetric_flow']*1e6, 
                  's-', linewidth=2, markersize=8, label=f'λ = {lambda_stoich}')
    ax2.set_xlabel('Current Density [A/cm²]')
    ax2.set_ylabel('Volumetric Flow Rate [cm³/s]')
    ax2.set_title('Flow Rate Requirements')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # 3. O₂ consumption vs current density
    ax3 = axes[1, 0]
    for lambda_stoich in stoichiometries:
        subset = op_df[op_df['stoichiometry'] == lambda_stoich]
        ax3.plot(subset['current_density'], subset['o2_consumption_rate']*1e6, 
                  '^-', linewidth=2, markersize=8, label=f'λ = {lambda_stoich}')
    ax3.set_xlabel('Current Density [A/cm²]')
    ax3.set_ylabel('O₂ Consumption Rate [mmol/s]')
    ax3.set_title('O₂ Consumption')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # 4. Bulk velocity vs current density
    ax4 = axes[1, 1]
    for lambda_stoich in stoichiometries:
        subset = op_df[op_df['stoichiometry'] == lambda_stoich]
        ax4.plot(subset['current_density'], subset['bulk_velocity'], 
                  'd-', linewidth=2, markersize=8, label=f'λ = {lambda_stoich}')
    ax4.set_xlabel('Current Density [A/cm²]')
    ax4.set_ylabel('Bulk Velocity [m/s]')
    ax4.set_title('Inlet Velocity')
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    
    plt.tight_layout()
    plt.show()
    
    # Save figure
    os.makedirs('../outputs/figures', exist_ok=True)
    plt.savefig('../outputs/figures/operating_conditions_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig('../outputs/figures/operating_conditions_analysis.pdf', bbox_inches='tight')
    print("\n✓ Analysis plots saved to outputs/figures/")
    
    # 8. Export Results
    print("\n8. EXPORTING RESULTS")
    print("-" * 30)
    
    # Export operating conditions to CSV
    output_dir = "../outputs/data"
    os.makedirs(output_dir, exist_ok=True)
    
    # Save operating conditions
    op_df.to_csv(f"{output_dir}/operating_conditions.csv", index=False)
    print(f"✓ Operating conditions exported to {output_dir}/operating_conditions.csv")
    
    # Create summary statistics
    summary_stats = {
        'total_cases': len(op_df),
        'current_density_range': f"{op_df['current_density'].min():.1f} - {op_df['current_density'].max():.1f} A/cm²",
        'stoichiometry_range': f"{op_df['stoichiometry'].min():.1f} - {op_df['stoichiometry'].max():.1f}",
        'reynolds_range': f"{op_df['reynolds_number'].min():.1f} - {op_df['reynolds_number'].max():.1f}",
        'flow_rate_range': f"{op_df['volumetric_flow'].min()*1e6:.2f} - {op_df['volumetric_flow'].max()*1e6:.2f} cm³/s"
    }
    
    print("\nStudy Summary:")
    for key, value in summary_stats.items():
        print(f"  {key.replace('_', ' ').title()}: {value}")
    
    # 9. Next Steps
    print("\n9. NEXT STEPS")
    print("-" * 30)
    print("This demo has shown the basic setup and analysis capabilities.")
    print("To run the full CFD simulation:")
    print("\n# Run single case")
    print("python src/main.py --mesh medium --current 1.2 --stoich 2.0")
    print("\n# Run mesh independence study")
    print("python src/main.py --study mesh_independence")
    print("\n# Run operating matrix study")
    print("python src/main.py --study operating_matrix")
    
    print("\nThe complete study will generate:")
    print("- Mesh files for all refinement levels")
    print("- CFD solutions for all operating conditions")
    print("- Publication-quality figures and data")
    print("- Comprehensive analysis reports")
    
    print("\n" + "="*60)
    print("DEMO COMPLETED SUCCESSFULLY!")
    print("="*60)

if __name__ == "__main__":
    main()
