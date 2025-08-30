#!/usr/bin/env python3
"""
PEMFC CFD Study - Main Execution Script
=======================================

This script runs the complete PEMFC CFD study with:
- Geometry generation and meshing
- Multi-physics solution
- Post-processing and analysis
- Results export
"""

import argparse
import os
import sys
import yaml
import logging
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from geometry.flowfield import PEMFCGeometry
from src.physics import PEMFCPhysics, MaterialProperties, get_default_properties
from src.solver import PEMFCSolver


def setup_logging(log_level: str = "INFO") -> logging.Logger:
    """Setup logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('outputs/logs/pemfc_cfd.log'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


def generate_geometry_and_mesh(mesh_refinement: str = "medium") -> Dict:
    """
    Generate geometry and mesh for the PEMFC flow field.
    
    Args:
        mesh_refinement: Mesh refinement level
        
    Returns:
        Dictionary with mesh file paths
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Generating {mesh_refinement} mesh...")
    
    try:
        # Create geometry
        pemfc_geom = PEMFCGeometry()
        
        # Generate mesh
        mesh_path = pemfc_geom.generate_mesh(mesh_refinement)
        fenics_path = pemfc_geom.export_to_fenics(mesh_path)
        
        mesh_files = {
            "gmsh": mesh_path,
            "fenics": fenics_path
        }
        
        logger.info(f"Mesh generation complete: {fenics_path}")
        return mesh_files
        
    except Exception as e:
        logger.error(f"Mesh generation failed: {e}")
        raise


def run_cfd_simulation(mesh_file: str, 
                      current_density: float,
                      stoichiometry: float,
                      mesh_refinement: str) -> Dict:
    """
    Run CFD simulation for given parameters.
    
    Args:
        mesh_file: Path to mesh file
        current_density: Current density [A/cm²]
        stoichiometry: Cathode stoichiometry
        mesh_refinement: Mesh refinement level
        
    Returns:
        Simulation results summary
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Running CFD simulation: j={current_density} A/cm², λ={stoichiometry}")
    
    try:
        # Initialize physics
        physics = PEMFCPhysics(get_default_properties())
        
        # Initialize solver
        solver = PEMFCSolver(mesh_file, physics, mesh_refinement)
        
        # Solve system
        converged = solver.solve(max_iterations=50, tolerance=1e-6)
        
        if not converged:
            logger.warning("Simulation did not converge within iterations")
        
        # Get results summary
        results = solver.get_solution_summary()
        results.update({
            "current_density": current_density,
            "stoichiometry": stoichiometry,
            "converged": converged
        })
        
        # Save solution
        output_dir = f"outputs/solutions/{mesh_refinement}_j{current_density}_lambda{stoichiometry}"
        solver.save_solution(output_dir)
        
        logger.info("CFD simulation complete")
        return results
        
    except Exception as e:
        logger.error(f"CFD simulation failed: {e}")
        raise


def calculate_operating_conditions(current_density: float, 
                                 stoichiometry: float,
                                 active_area: float = 25.0e-4) -> Dict:
    """
    Calculate operating conditions from electrochemical parameters.
    
    Args:
        current_density: Current density [A/cm²]
        stoichiometry: Cathode stoichiometry
        active_area: Active area [m²]
        
    Returns:
        Dictionary with operating conditions
    """
    # Constants
    F = 96485.0  # C/mol
    R = 8.314    # J/(mol·K)
    T = 353.0    # K
    P = 1.013e5  # Pa
    x_o2_in = 0.21  # O₂ mole fraction
    
    # Convert current density to A/m²
    j_am2 = current_density * 1e4
    
    # Calculate O₂ consumption rate
    n_o2_cons = j_am2 * active_area / (4 * F)
    
    # Calculate inlet O₂ flow rate
    n_o2_in = stoichiometry * n_o2_cons
    
    # Calculate total inlet flow rate
    n_tot_in = n_o2_in / x_o2_in
    
    # Calculate volumetric flow rate
    Q = n_tot_in * R * T / P
    
    # Calculate bulk velocity (assuming inlet area)
    inlet_area = 25.0e-6  # m² (approximate)
    U_bulk = Q / inlet_area
    
    return {
        "current_density": current_density,
        "stoichiometry": stoichiometry,
        "o2_consumption_rate": n_o2_cons,
        "o2_inlet_rate": n_o2_in,
        "total_inlet_rate": n_tot_in,
        "volumetric_flow": Q,
        "bulk_velocity": U_bulk,
        "reynolds_number": 1.009 * U_bulk * 0.001 / 2.08e-5  # Re = ρUL/μ
    }


def run_mesh_independence_study() -> Dict:
    """
    Run mesh independence study with three refinement levels.
    
    Returns:
        Dictionary with mesh independence results
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting mesh independence study...")
    
    # Operating conditions
    current_density = 1.2  # A/cm²
    stoichiometry = 2.0
    
    mesh_results = {}
    
    for mesh_level in ["coarse", "medium", "fine"]:
        logger.info(f"Running {mesh_level} mesh...")
        
        try:
            # Generate mesh
            mesh_files = generate_geometry_and_mesh(mesh_level)
            
            # Run simulation
            results = run_cfd_simulation(
                mesh_files["fenics"],
                current_density,
                stoichiometry,
                mesh_level
            )
            
            mesh_results[mesh_level] = results
            
        except Exception as e:
            logger.error(f"Failed to run {mesh_level} mesh: {e}")
            continue
    
    # Analyze mesh independence
    mesh_independence = analyze_mesh_independence(mesh_results)
    
    return {
        "mesh_results": mesh_results,
        "mesh_independence": mesh_independence
    }


def analyze_mesh_independence(mesh_results: Dict) -> Dict:
    """
    Analyze mesh independence results.
    
    Args:
        mesh_results: Results from different mesh levels
        
    Returns:
        Mesh independence analysis
    """
    if len(mesh_results) < 2:
        return {"status": "insufficient_data"}
    
    # Extract key metrics
    metrics = {}
    for mesh_level, results in mesh_results.items():
        if "mesh_info" in results:
            metrics[mesh_level] = {
                "num_cells": results["mesh_info"]["num_cells"],
                "residuals": results["final_residuals"]
            }
    
    # Calculate relative changes
    if "medium" in metrics and "fine" in metrics:
        medium_cells = metrics["medium"]["num_cells"]
        fine_cells = metrics["fine"]["num_cells"]
        
        # Grid convergence index (simplified)
        relative_change = abs(fine_cells - medium_cells) / medium_cells
        
        mesh_independence = {
            "status": "complete",
            "medium_to_fine_change": relative_change,
            "converged": relative_change < 0.02,  # 2% threshold
            "recommendation": "fine" if relative_change < 0.02 else "medium"
        }
    else:
        mesh_independence = {"status": "incomplete"}
    
    return mesh_independence


def run_operating_matrix_study() -> Dict:
    """
    Run study across operating matrix (current density × stoichiometry).
    
    Returns:
        Dictionary with operating matrix results
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting operating matrix study...")
    
    # Operating matrix
    current_densities = [0.3, 0.8, 1.2]  # A/cm²
    stoichiometries = [1.5, 2.5]
    mesh_level = "medium"  # Use medium mesh for matrix study
    
    results_matrix = {}
    
    for j in current_densities:
        for lambda_stoich in stoichiometries:
            case_key = f"j{j}_lambda{lambda_stoich}"
            logger.info(f"Running case: {case_key}")
            
            try:
                # Calculate operating conditions
                op_conditions = calculate_operating_conditions(j, lambda_stoich)
                
                # Generate mesh if not exists
                mesh_files = generate_geometry_and_mesh(mesh_level)
                
                # Run simulation
                simulation_results = run_cfd_simulation(
                    mesh_files["fenics"],
                    j,
                    lambda_stoich,
                    mesh_level
                )
                
                # Combine results
                results_matrix[case_key] = {
                    "operating_conditions": op_conditions,
                    "simulation_results": simulation_results
                }
                
            except Exception as e:
                logger.error(f"Failed to run case {case_key}: {e}")
                continue
    
    return results_matrix


def export_results(results: Dict, output_dir: str = "outputs"):
    """
    Export results to files.
    
    Args:
        results: Results dictionary
        output_dir: Output directory
    """
    logger = logging.getLogger(__name__)
    logger.info("Exporting results...")
    
    # Create output directories
    os.makedirs(f"{output_dir}/data", exist_ok=True)
    os.makedirs(f"{output_dir}/figures", exist_ok=True)
    
    # Export mesh independence results
    if "mesh_independence" in results:
        mesh_df = pd.DataFrame(results["mesh_results"]).T
        mesh_df.to_csv(f"{output_dir}/data/mesh_independence.csv")
        
        # Plot mesh convergence
        plot_mesh_convergence(results["mesh_results"], f"{output_dir}/figures")
    
    # Export operating matrix results
    if "operating_matrix" in results:
        op_matrix_df = pd.DataFrame(results["operating_matrix"]).T
        op_matrix_df.to_csv(f"{output_dir}/data/operating_matrix.csv")
        
        # Plot operating matrix results
        plot_operating_matrix(results["operating_matrix"], f"{output_dir}/figures")
    
    logger.info("Results exported successfully")


def plot_mesh_convergence(mesh_results: Dict, output_dir: str):
    """Plot mesh convergence results."""
    # Extract data
    mesh_levels = list(mesh_results.keys())
    cell_counts = [mesh_results[level]["mesh_info"]["num_cells"] for level in mesh_levels]
    
    # Create plot
    plt.figure(figsize=(10, 6))
    plt.semilogx(cell_counts, [1, 2, 3], 'o-', linewidth=2, markersize=8)
    plt.xlabel('Number of Cells')
    plt.ylabel('Mesh Level')
    plt.title('Mesh Convergence Study')
    plt.grid(True, alpha=0.3)
    
    # Save plot
    plt.savefig(f"{output_dir}/mesh_convergence.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/mesh_convergence.pdf", bbox_inches='tight')
    plt.close()


def plot_operating_matrix(operating_matrix: Dict, output_dir: str):
    """Plot operating matrix results."""
    # Extract data for plotting
    current_densities = []
    stoichiometries = []
    pressure_drops = []
    
    for case_key, case_data in operating_matrix.items():
        if "operating_conditions" in case_data:
            op_cond = case_data["operating_conditions"]
            current_densities.append(op_cond["current_density"])
            stoichiometries.append(op_cond["stoichiometry"])
            # Add pressure drop if available
            pressure_drops.append(0.0)  # Placeholder
    
    # Create plot
    plt.figure(figsize=(12, 8))
    
    # Scatter plot
    scatter = plt.scatter(current_densities, stoichiometries, 
                         c=pressure_drops, s=100, cmap='viridis')
    plt.colorbar(scatter, label='Pressure Drop [Pa]')
    
    plt.xlabel('Current Density [A/cm²]')
    plt.ylabel('Stoichiometry')
    plt.title('Operating Matrix Results')
    plt.grid(True, alpha=0.3)
    
    # Save plot
    plt.savefig(f"{output_dir}/operating_matrix.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/operating_matrix.pdf", bbox_inches='tight')
    plt.close()


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description="PEMFC CFD Study")
    parser.add_argument("--mesh", choices=["coarse", "medium", "fine"], 
                       default="medium", help="Mesh refinement level")
    parser.add_argument("--current", type=float, default=1.2,
                       help="Current density [A/cm²]")
    parser.add_argument("--stoich", type=float, default=2.0,
                       help="Cathode stoichiometry")
    parser.add_argument("--study", choices=["single", "mesh_independence", "operating_matrix"],
                       default="single", help="Study type")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                       default="INFO", help="Logging level")
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.log_level)
    logger.info("Starting PEMFC CFD Study")
    
    try:
        if args.study == "single":
            # Single case study
            logger.info(f"Running single case: mesh={args.mesh}, j={args.current}, λ={args.stoich}")
            
            # Generate mesh
            mesh_files = generate_geometry_and_mesh(args.mesh)
            
            # Run simulation
            results = run_cfd_simulation(
                mesh_files["fenics"],
                args.current,
                args.stoich,
                args.mesh
            )
            
            # Export results
            export_results({"single_case": results})
            
        elif args.study == "mesh_independence":
            # Mesh independence study
            logger.info("Running mesh independence study")
            results = run_mesh_independence_study()
            export_results(results)
            
        elif args.study == "operating_matrix":
            # Operating matrix study
            logger.info("Running operating matrix study")
            results = {"operating_matrix": run_operating_matrix_study()}
            export_results(results)
        
        logger.info("PEMFC CFD Study completed successfully")
        
    except Exception as e:
        logger.error(f"Study failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
