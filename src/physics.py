"""
PEMFC Physics Module
====================

This module contains the governing equations and constitutive laws for the PEMFC CFD study:
- Incompressible Navier-Stokes with Brinkman-Forchheimer porous media
- Energy equation (convection + conduction)
- Multi-species transport (O₂, H₂O)
- Butler-Volmer electrochemistry
"""

import numpy as np
from typing import Dict, Tuple, Optional
from dataclasses import dataclass


@dataclass
class MaterialProperties:
    """Material properties for PEMFC components."""
    
    # Fluid properties (air at 353K, 1atm)
    rho_air: float = 1.009  # kg/m³
    mu_air: float = 2.08e-5  # Pa·s
    cp_air: float = 1008.0  # J/(kg·K)
    k_air: float = 0.0301  # W/(m·K)
    
    # GDL properties
    gdl_porosity: float = 0.7
    gdl_permeability: float = 1.0e-12  # m²
    gdl_forchheimer_coeff: float = 0.1
    gdl_tortuosity: float = 1.5
    k_gdl_axial: float = 20.0  # W/(m·K)
    k_gdl_transverse: float = 1.0  # W/(m·K)
    
    # Membrane properties
    k_membrane: float = 0.95  # W/(m·K)
    membrane_conductivity: float = 10.0  # S/m
    
    # Electrochemical constants
    faraday_constant: float = 96485.0  # C/mol
    gas_constant: float = 8.314  # J/(mol·K)
    alpha_anodic: float = 0.5
    alpha_cathodic: float = 0.5
    exchange_current_density: float = 1.0e-3  # A/m²
    reversible_voltage: float = 1.23  # V
    
    # Species properties
    d_o2_air: float = 2.0e-5  # m²/s
    d_h2o_air: float = 2.8e-5  # m²/s


class PEMFCPhysics:
    """PEMFC physics solver with all governing equations."""
    
    def __init__(self, props: MaterialProperties):
        """
        Initialize physics solver.
        
        Args:
            props: Material properties
        """
        self.props = props
        
    def brinkman_forchheimer_coeffs(self, domain: str) -> Tuple[float, float]:
        """
        Get Brinkman-Forchheimer coefficients for porous media.
        
        Args:
            domain: Domain identifier ("FLUID", "GDL", "MEMBRANE")
            
        Returns:
            Tuple of (alpha, beta) coefficients
        """
        if domain == "FLUID":
            return 0.0, 0.0
        elif domain == "GDL":
            alpha = self.props.mu_air / self.props.gdl_permeability
            beta = self.props.rho_air * self.props.gdl_forchheimer_coeff / 2.0
            return alpha, beta
        elif domain == "MEMBRANE":
            # Membrane is solid, no flow
            return 1e20, 0.0
        else:
            raise ValueError(f"Unknown domain: {domain}")
    
    def effective_diffusivity(self, species: str, domain: str) -> float:
        """
        Calculate effective diffusivity in porous media.
        
        Args:
            species: Species identifier ("O2", "H2O")
            domain: Domain identifier
            
        Returns:
            Effective diffusivity [m²/s]
        """
        if domain == "FLUID":
            if species == "O2":
                return self.props.d_o2_air
            elif species == "H2O":
                return self.props.d_h2o_air
            else:
                raise ValueError(f"Unknown species: {species}")
        
        elif domain in ["GDL", "MEMBRANE"]:
            # Bruggeman correlation
            if species == "O2":
                base_diff = self.props.d_o2_air
            elif species == "H2O":
                base_diff = self.props.d_h2o_air
            else:
                raise ValueError(f"Unknown species: {species}")
            
            if domain == "GDL":
                porosity = self.props.gdl_porosity
                tortuosity = self.props.gdl_tortuosity
            else:  # MEMBRANE
                porosity = 0.1  # Nafion water content
                tortuosity = 2.0
            
            return base_diff * (porosity ** tortuosity)
        
        else:
            raise ValueError(f"Unknown domain: {domain}")
    
    def butler_volmer_current(self, 
                             eta: float, 
                             c_o2: float, 
                             c_o2_ref: float,
                             temperature: float) -> float:
        """
        Calculate Butler-Volmer current density.
        
        Args:
            eta: Overpotential [V]
            c_o2: Local O₂ concentration [mol/m³]
            c_o2_ref: Reference O₂ concentration [mol/m³]
            temperature: Local temperature [K]
            
        Returns:
            Current density [A/m²]
        """
        # Temperature-dependent exchange current density
        t_ref = 353.0  # K
        j0_temp = self.props.exchange_current_density * np.exp(
            -self.props.gas_constant * (1/temperature - 1/t_ref) / 
            (self.props.gas_constant * t_ref)
        )
        
        # Concentration dependence
        j0 = j0_temp * (c_o2 / c_o2_ref) ** 0.5
        
        # Butler-Volmer equation
        j = j0 * (
            np.exp(self.props.alpha_anodic * self.props.faraday_constant * eta / 
                   (self.props.gas_constant * temperature)) -
            np.exp(-self.props.alpha_cathodic * self.props.faraday_constant * eta / 
                   (self.props.gas_constant * temperature))
        )
        
        return j
    
    def electrochemical_sources(self, 
                               j: float, 
                               domain: str) -> Dict[str, float]:
        """
        Calculate electrochemical source terms.
        
        Args:
            j: Current density [A/m²]
            domain: Domain identifier
            
        Returns:
            Dictionary of source terms for each species and energy
        """
        sources = {
            "O2": 0.0,
            "H2O": 0.0,
            "heat_reaction": 0.0,
            "heat_ohmic": 0.0
        }
        
        if domain == "GDL":
            # O₂ consumption
            sources["O2"] = -j / (4.0 * self.props.faraday_constant)
            
            # H₂O production
            sources["H2O"] = j / (2.0 * self.props.faraday_constant)
            
            # Reaction heat
            sources["heat_reaction"] = j * self.props.reversible_voltage
            
            # Ohmic heat (simplified)
            sources["heat_ohmic"] = j**2 / self.props.membrane_conductivity
        
        return sources
    
    def calculate_reynolds_number(self, velocity: float, length: float) -> float:
        """
        Calculate Reynolds number.
        
        Args:
            velocity: Characteristic velocity [m/s]
            length: Characteristic length [m]
            
        Returns:
            Reynolds number
        """
        return self.props.rho_air * velocity * length / self.props.mu_air
    
    def calculate_peclet_number(self, velocity: float, length: float, 
                               diffusivity: float) -> float:
        """
        Calculate Peclet number.
        
        Args:
            velocity: Characteristic velocity [m/s]
            length: Characteristic length [m]
            diffusivity: Diffusivity [m²/s]
            
        Returns:
            Peclet number
        """
        return velocity * length / diffusivity
    
    def calculate_pumping_power(self, pressure_drop: float, 
                               volumetric_flow: float) -> float:
        """
        Calculate pumping power.
        
        Args:
            pressure_drop: Pressure drop [Pa]
            volumetric_flow: Volumetric flow rate [m³/s]
            
        Returns:
            Pumping power [W]
        """
        return pressure_drop * volumetric_flow


def get_default_properties() -> MaterialProperties:
    """Get default material properties for PEMFC."""
    return MaterialProperties()


def validate_physics_inputs(velocity: float, temperature: float, 
                           pressure: float, o2_mole_fraction: float) -> bool:
    """
    Validate physics input parameters.
    
    Args:
        velocity: Inlet velocity [m/s]
        temperature: Temperature [K]
        pressure: Pressure [Pa]
        o2_mole_fraction: O₂ mole fraction
        
    Returns:
        True if inputs are valid
    """
    if velocity <= 0 or velocity > 100:
        return False
    
    if temperature < 273 or temperature > 373:
        return False
    
    if pressure < 0.5e5 or pressure > 2e5:
        return False
    
    if o2_mole_fraction <= 0 or o2_mole_fraction > 1:
        return False
    
    return True
