"""
PEMFC FEniCSx Solver
=====================

This module implements the main FEniCSx solver for the PEMFC CFD study with:
- Incompressible Navier-Stokes with porous media
- Energy equation
- Multi-species transport
- Electrochemical coupling
"""

import numpy as np
import dolfinx
from dolfinx import fem, mesh, io
from dolfinx.fem import Function, FunctionSpace, Constant
from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from petsc4py import PETSc
import ufl
from typing import Dict, List, Tuple, Optional
import logging
import os

from .physics import PEMFCPhysics, MaterialProperties


class PEMFCSolver:
    """Main PEMFC CFD solver using FEniCSx."""
    
    def __init__(self, 
                 mesh_file: str,
                 physics: PEMFCPhysics,
                 mesh_refinement: str = "medium"):
        """
        Initialize PEMFC solver.
        
        Args:
            mesh_file: Path to mesh file (.xdmf)
            physics: Physics module instance
            mesh_refinement: Mesh refinement level
        """
        self.mesh_file = mesh_file
        self.physics = physics
        self.mesh_refinement = mesh_refinement
        
        # Load mesh
        self.mesh = self._load_mesh()
        
        # Initialize function spaces
        self._setup_function_spaces()
        
        # Initialize functions
        self._setup_functions()
        
        # Setup boundary conditions
        self._setup_boundary_conditions()
        
        # Setup variational forms
        self._setup_variational_forms()
        
        # Setup solvers
        self._setup_solvers()
        
        # Convergence tracking
        self.residuals = {
            "velocity": [],
            "pressure": [],
            "temperature": [],
            "o2": [],
            "h2o": []
        }
        
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
    
    def _load_mesh(self) -> dolfinx.mesh.Mesh:
        """Load mesh from file."""
        try:
            with dolfinx.io.XDMFFile(self.mesh.comm, self.mesh_file, "r") as xdmf:
                mesh = xdmf.read_mesh()
            self.logger.info(f"Loaded mesh: {self.mesh_file}")
            return mesh
        except Exception as e:
            self.logger.error(f"Failed to load mesh: {e}")
            raise
    
    def _setup_function_spaces(self):
        """Setup function spaces for all variables."""
        # Velocity-pressure space (Taylor-Hood P2/P1)
        self.V = dolfinx.fem.VectorFunctionSpace(self.mesh, ("CG", 2))
        self.Q = dolfinx.fem.FunctionSpace(self.mesh, ("CG", 1))
        self.W = dolfinx.fem.FunctionSpace(self.mesh, ("CG", 2), dim=3)
        
        # Scalar spaces
        self.T_space = dolfinx.fem.FunctionSpace(self.mesh, ("CG", 2))
        self.O2_space = dolfinx.fem.FunctionSpace(self.mesh, ("CG", 2))
        self.H2O_space = dolfinx.fem.FunctionSpace(self.mesh, ("CG", 2))
        
        # Mixed space for velocity-pressure coupling
        self.VQ = self.V * self.Q
        
        self.logger.info("Function spaces initialized")
    
    def _setup_functions(self):
        """Initialize solution functions."""
        # Current solution
        self.u = Function(self.V)  # Velocity
        self.p = Function(self.Q)  # Pressure
        self.T = Function(self.T_space)  # Temperature
        self.Y_O2 = Function(self.O2_space)  # O₂ mass fraction
        self.Y_H2O = Function(self.H2O_space)  # H₂O mass fraction
        
        # Previous solution (for time stepping)
        self.u_prev = Function(self.V)
        self.T_prev = Function(self.T_space)
        self.Y_O2_prev = Function(self.O2_space)
        self.Y_H2O_prev = Function(self.H2O_space)
        
        # Test functions
        self.v = ufl.TestFunction(self.V)
        self.q = ufl.TestFunction(self.Q)
        self.w_T = ufl.TestFunction(self.T_space)
        self.w_O2 = ufl.TestFunction(self.O2_space)
        self.w_H2O = ufl.TestFunction(self.H2O_space)
        
        # Trial functions
        self.u_trial = ufl.TrialFunction(self.V)
        self.p_trial = ufl.TrialFunction(self.Q)
        self.T_trial = ufl.TrialFunction(self.T_space)
        self.Y_O2_trial = ufl.TrialFunction(self.O2_space)
        self.Y_H2O_trial = ufl.TrialFunction(self.H2O_space)
        
        self.logger.info("Solution functions initialized")
    
    def _setup_boundary_conditions(self):
        """Setup boundary conditions."""
        # Get boundary facets
        boundary_facets = dolfinx.mesh.locate_entities_boundary(
            self.mesh, self.mesh.topology.dim - 1,
            lambda x: np.logical_or(x[0] < 1e-10, x[0] > 1e-10)
        )
        
        # Inlet boundary (x = 0)
        inlet_facets = dolfinx.mesh.locate_entities_boundary(
            self.mesh, self.mesh.topology.dim - 1,
            lambda x: x[0] < 1e-10
        )
        
        # Outlet boundary (x = L)
        outlet_facets = dolfinx.mesh.locate_entities_boundary(
            self.mesh, self.mesh.topology.dim - 1,
            lambda x: x[0] > 1e-10
        )
        
        # Boundary conditions
        self.bcs = {
            "velocity": [],
            "pressure": [],
            "temperature": [],
            "o2": [],
            "h2o": []
        }
        
        # Velocity BCs
        inlet_velocity = Constant(self.mesh, PETSc.ScalarType([1.0, 0.0, 0.0]))
        inlet_bc = dolfinx.fem.dirichletbc(inlet_velocity, inlet_facets, self.V)
        self.bcs["velocity"].append(inlet_bc)
        
        # Pressure BC (outlet reference)
        outlet_pressure = Constant(self.mesh, PETSc.ScalarType(0.0))
        outlet_bc = dolfinx.fem.dirichletbc(outlet_pressure, outlet_facets, self.Q)
        self.bcs["pressure"].append(outlet_bc)
        
        # Temperature BCs
        inlet_temp = Constant(self.mesh, PETSc.ScalarType(353.0))
        inlet_temp_bc = dolfinx.fem.dirichletbc(inlet_temp, inlet_facets, self.T_space)
        self.bcs["temperature"].append(inlet_temp_bc)
        
        # Species BCs
        inlet_o2 = Constant(self.mesh, PETSc.ScalarType(0.21))
        inlet_h2o = Constant(self.mesh, PETSc.ScalarType(0.0))
        inlet_o2_bc = dolfinx.fem.dirichletbc(inlet_o2, inlet_facets, self.O2_space)
        inlet_h2o_bc = dolfinx.fem.dirichletbc(inlet_h2o, inlet_facets, self.H2O_space)
        self.bcs["o2"].append(inlet_o2_bc)
        self.bcs["h2o"].append(inlet_h2o_bc)
        
        self.logger.info("Boundary conditions initialized")
    
    def _setup_variational_forms(self):
        """Setup variational forms for all equations."""
        # Material properties
        rho = Constant(self.mesh, PETSc.ScalarType(self.physics.props.rho_air))
        mu = Constant(self.mesh, PETSc.ScalarType(self.physics.props.mu_air))
        cp = Constant(self.mesh, PETSc.ScalarType(self.physics.props.cp_air))
        k = Constant(self.mesh, PETSc.ScalarType(self.physics.props.k_air))
        
        # Navier-Stokes form
        self._setup_navier_stokes_form(rho, mu)
        
        # Energy form
        self._setup_energy_form(rho, cp, k)
        
        # Species transport forms
        self._setup_species_forms(rho)
        
        self.logger.info("Variational forms initialized")
    
    def _setup_navier_stokes_form(self, rho: Constant, mu: Constant):
        """Setup Navier-Stokes variational form."""
        # Convection term
        conv = rho * ufl.dot(ufl.dot(self.u, ufl.grad(self.u)), self.v)
        
        # Viscous term
        viscous = 2 * mu * ufl.inner(ufl.grad(self.u), ufl.grad(self.v))
        
        # Pressure term
        pressure = -self.p * ufl.div(self.v)
        
        # Continuity
        continuity = self.q * ufl.div(self.u)
        
        # Brinkman-Forchheimer terms (will be added per domain)
        self.ns_form = conv + viscous + pressure + continuity
        
        # LHS for linear system
        self.ns_lhs = ufl.inner(ufl.grad(self.u_trial), ufl.grad(self.v)) + \
                      ufl.inner(self.p_trial, ufl.div(self.v)) + \
                      ufl.inner(ufl.div(self.u_trial), self.q)
    
    def _setup_energy_form(self, rho: Constant, cp: Constant, k: Constant):
        """Setup energy equation variational form."""
        # Convection
        conv = rho * cp * ufl.dot(self.u, ufl.grad(self.T)) * self.w_T
        
        # Conduction
        conduction = k * ufl.dot(ufl.grad(self.T), ufl.grad(self.w_T))
        
        # Source terms (will be added)
        self.energy_form = conv + conduction
        
        # LHS
        self.energy_lhs = ufl.inner(ufl.grad(self.T_trial), ufl.grad(self.w_T))
    
    def _setup_species_forms(self, rho: Constant):
        """Setup species transport variational forms."""
        # O₂ transport
        conv_o2 = rho * ufl.dot(self.u, ufl.grad(self.Y_O2)) * self.w_O2
        diff_o2 = rho * self.physics.props.d_o2_air * ufl.dot(ufl.grad(self.Y_O2), ufl.grad(self.w_O2))
        self.o2_form = conv_o2 + diff_o2
        
        # H₂O transport
        conv_h2o = rho * ufl.dot(self.u, ufl.grad(self.Y_H2O)) * self.w_H2O
        diff_h2o = rho * self.physics.props.d_h2o_air * ufl.dot(ufl.grad(self.Y_H2O), ufl.grad(self.w_H2O))
        self.h2o_form = conv_h2o + diff_h2o
        
        # LHS forms
        self.o2_lhs = ufl.inner(ufl.grad(self.Y_O2_trial), ufl.grad(self.w_O2))
        self.h2o_lhs = ufl.inner(ufl.grad(self.Y_H2O_trial), ufl.grad(self.w_H2O))
    
    def _setup_solvers(self):
        """Setup PETSc solvers."""
        # Navier-Stokes solver
        self.ns_problem = LinearProblem(
            self.ns_lhs, 
            self.ns_form, 
            bcs=self.bcs["velocity"] + self.bcs["pressure"],
            petsc_options={"ksp_type": "gmres", "pc_type": "fieldsplit"}
        )
        
        # Energy solver
        self.energy_problem = LinearProblem(
            self.energy_lhs,
            self.energy_form,
            bcs=self.bcs["temperature"],
            petsc_options={"ksp_type": "cg", "pc_type": "hypre"}
        )
        
        # Species solvers
        self.o2_problem = LinearProblem(
            self.o2_lhs,
            self.o2_form,
            bcs=self.bcs["o2"],
            petsc_options={"ksp_type": "cg", "pc_type": "hypre"}
        )
        
        self.h2o_problem = LinearProblem(
            self.h2o_lhs,
            self.h2o_form,
            bcs=self.bcs["h2o"],
            petsc_options={"ksp_type": "cg", "pc_type": "hypre"}
        )
        
        self.logger.info("PETSc solvers initialized")
    
    def solve(self, max_iterations: int = 50, tolerance: float = 1e-6) -> bool:
        """
        Solve the coupled system.
        
        Args:
            max_iterations: Maximum Picard iterations
            tolerance: Convergence tolerance
            
        Returns:
            True if converged
        """
        self.logger.info("Starting coupled solution...")
        
        for iteration in range(max_iterations):
            self.logger.info(f"Iteration {iteration + 1}")
            
            # Solve Navier-Stokes
            converged_ns = self._solve_navier_stokes()
            if not converged_ns:
                self.logger.error("Navier-Stokes failed to converge")
                return False
            
            # Solve energy equation
            converged_energy = self._solve_energy()
            if not converged_energy:
                self.logger.error("Energy equation failed to converge")
                return False
            
            # Solve species transport
            converged_o2 = self._solve_o2()
            converged_h2o = self._solve_h2o()
            if not (converged_o2 and converged_h2o):
                self.logger.error("Species transport failed to converge")
                return False
            
            # Check convergence
            if self._check_convergence(tolerance):
                self.logger.info(f"Converged in {iteration + 1} iterations")
                return True
            
            # Update previous solutions
            self._update_previous_solutions()
        
        self.logger.warning("Failed to converge within maximum iterations")
        return False
    
    def _solve_navier_stokes(self) -> bool:
        """Solve Navier-Stokes equations."""
        try:
            solution = self.ns_problem.solve()
            self.u.x.array[:] = solution.x.array[:self.V.dofmap.index_map.size_local * self.V.dofmap.index_map_bs]
            self.p.x.array[:] = solution.x.array[self.V.dofmap.index_map.size_local * self.V.dofmap.index_map_bs:]
            return True
        except Exception as e:
            self.logger.error(f"Navier-Stokes solve failed: {e}")
            return False
    
    def _solve_energy(self) -> bool:
        """Solve energy equation."""
        try:
            solution = self.energy_problem.solve()
            self.T.x.array[:] = solution.x.array[:]
            return True
        except Exception as e:
            self.logger.error(f"Energy solve failed: {e}")
            return False
    
    def _solve_o2(self) -> bool:
        """Solve O₂ transport equation."""
        try:
            solution = self.o2_problem.solve()
            self.Y_O2.x.array[:] = solution.x.array[:]
            return True
        except Exception as e:
            self.logger.error(f"O₂ solve failed: {e}")
            return False
    
    def _solve_h2o(self) -> bool:
        """Solve H₂O transport equation."""
        try:
            solution = self.h2o_problem.solve()
            self.Y_H2O.x.array[:] = solution.x.array[:]
            return True
        except Exception as e:
            self.logger.error(f"H₂O solve failed: {e}")
            return False
    
    def _check_convergence(self, tolerance: float) -> bool:
        """Check if solution has converged."""
        # Calculate residuals
        residual_u = self._calculate_velocity_residual()
        residual_p = self._calculate_pressure_residual()
        residual_T = self._calculate_temperature_residual()
        residual_o2 = self._calculate_o2_residual()
        residual_h2o = self._calculate_h2o_residual()
        
        # Store residuals
        self.residuals["velocity"].append(residual_u)
        self.residuals["pressure"].append(residual_p)
        self.residuals["temperature"].append(residual_T)
        self.residuals["o2"].append(residual_o2)
        self.residuals["h2o"].append(residual_h2o)
        
        # Check convergence
        max_residual = max(residual_u, residual_p, residual_T, residual_o2, residual_h2o)
        converged = max_residual < tolerance
        
        self.logger.info(f"Max residual: {max_residual:.2e}")
        
        return converged
    
    def _calculate_velocity_residual(self) -> float:
        """Calculate velocity residual."""
        # Simplified residual calculation
        return np.linalg.norm(self.u.x.array - self.u_prev.x.array)
    
    def _calculate_pressure_residual(self) -> float:
        """Calculate pressure residual."""
        return np.linalg.norm(self.p.x.array)
    
    def _calculate_temperature_residual(self) -> float:
        """Calculate temperature residual."""
        return np.linalg.norm(self.T.x.array - self.T_prev.x.array)
    
    def _calculate_o2_residual(self) -> float:
        """Calculate O₂ residual."""
        return np.linalg.norm(self.Y_O2.x.array - self.Y_O2_prev.x.array)
    
    def _calculate_h2o_residual(self) -> float:
        """Calculate H₂O residual."""
        return np.linalg.norm(self.Y_H2O.x.array - self.Y_H2O_prev.x.array)
    
    def _update_previous_solutions(self):
        """Update previous solution values."""
        self.u_prev.x.array[:] = self.u.x.array[:]
        self.T_prev.x.array[:] = self.T.x.array[:]
        self.Y_O2_prev.x.array[:] = self.Y_O2.x.array[:]
        self.Y_H2O_prev.x.array[:] = self.Y_H2O.x.array[:]
    
    def save_solution(self, output_dir: str):
        """Save solution to files."""
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Save mesh
        with dolfinx.io.VTXWriter(self.mesh.comm, f"{output_dir}/mesh.bp", [self.mesh]) as vtx:
            vtx.write(0.0)
        
        # Save solution fields
        with dolfinx.io.VTXWriter(self.mesh.comm, f"{output_dir}/solution.bp", 
                                 [self.u, self.p, self.T, self.Y_O2, self.Y_H2O]) as vtx:
            vtx.write(0.0)
        
        self.logger.info(f"Solution saved to {output_dir}")
    
    def get_solution_summary(self) -> Dict:
        """Get summary of solution results."""
        return {
            "mesh_refinement": self.mesh_refinement,
            "final_residuals": {
                "velocity": self.residuals["velocity"][-1] if self.residuals["velocity"] else None,
                "pressure": self.residuals["pressure"][-1] if self.residuals["pressure"] else None,
                "temperature": self.residuals["temperature"][-1] if self.residuals["temperature"] else None,
                "o2": self.residuals["o2"][-1] if self.residuals["o2"] else None,
                "h2o": self.residuals["h2o"][-1] if self.residuals["h2o"] else None
            },
            "convergence_history": self.residuals,
            "mesh_info": {
                "num_cells": self.mesh.topology.index_map(self.mesh.topology.dim).size_global,
                "num_vertices": self.mesh.topology.index_map(0).size_global
            }
        }
