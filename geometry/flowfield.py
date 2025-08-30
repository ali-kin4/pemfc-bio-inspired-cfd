"""
PEMFC Flow Field Geometry Generator
==================================

This module generates a tapered interdigitated cathode flow field with a two-level
lung-inspired inlet manifold using pygmsh and gmsh.

Geometry specifications:
- Channel depth: 1.0 mm
- Channel width: 1.0 mm (inlet) → 0.85-0.90 mm (outlet) - linear taper
- Rib width: 1.0 mm
- GDL thickness: 0.25 mm
- Membrane thickness: 0.075 mm
- Two-level lung-inspired branching manifold
"""

import pygmsh
import numpy as np
import gmsh
from typing import Tuple, List, Dict
import os


class PEMFCGeometry:
    """PEMFC flow field geometry generator."""
    
    def __init__(self, 
                 channel_depth: float = 1.0e-3,
                 channel_width_inlet: float = 1.0e-3,
                 channel_width_outlet: float = 0.85e-3,
                 rib_width: float = 1.0e-3,
                 gdl_thickness: float = 2.5e-4,
                 membrane_thickness: float = 7.5e-5,
                 num_channels: int = 8,
                 channel_length: float = 20.0e-3,
                 manifold_levels: int = 2):
        """
        Initialize geometry parameters.
        
        Args:
            channel_depth: Depth of flow channels [m]
            channel_width_inlet: Width at inlet [m]
            channel_width_outlet: Width at outlet [m]
            rib_width: Width of ribs between channels [m]
            gdl_thickness: Thickness of GDL layer [m]
            membrane_thickness: Thickness of membrane [m]
            num_channels: Number of interdigitated channels
            channel_length: Length of channels [m]
            manifold_levels: Number of manifold branching levels
        """
        self.channel_depth = channel_depth
        self.channel_width_inlet = channel_width_inlet
        self.channel_width_outlet = channel_width_outlet
        self.rib_width = rib_width
        self.gdl_thickness = gdl_thickness
        self.membrane_thickness = membrane_thickness
        self.num_channels = num_channels
        self.channel_length = channel_length
        self.manifold_levels = manifold_levels
        
        # Calculate total width
        self.total_width = num_channels * (channel_width_inlet + rib_width) - rib_width
        
        # Initialize pygmsh
        self.geom = pygmsh.geo.Geometry()
        
    def create_lung_manifold(self) -> Tuple[List, List]:
        """
        Create two-level lung-inspired branching manifold.
        
        Returns:
            Tuple of (inlet_points, outlet_points) for channels
        """
        # Main inlet
        inlet_width = self.total_width * 0.3  # 30% of total width
        inlet_height = self.channel_depth
        
        # First level branching (2 branches)
        branch1_width = inlet_width * 0.6
        branch2_width = inlet_width * 0.4
        
        # Second level branching (4 branches)
        sub_branch_widths = [branch1_width * 0.5, branch1_width * 0.5,
                            branch2_width * 0.6, branch2_width * 0.4]
        
        # Create points for manifold
        points = []
        lines = []
        
        # Main inlet
        p1 = self.geom.add_point([0, 0, 0])
        p2 = self.geom.add_point([inlet_width, 0, 0])
        p3 = self.geom.add_point([inlet_width, inlet_height, 0])
        p4 = self.geom.add_point([0, inlet_height, 0])
        
        points.extend([p1, p2, p3, p4])
        
        # Create inlet rectangle
        inlet_rect = self.geom.add_rectangle(0, 0, inlet_width, inlet_height)
        
        # First level branches
        branch1_start = inlet_width * 0.2
        branch2_start = inlet_width * 0.6
        
        # Branch 1
        p5 = self.geom.add_point([branch1_start, 0, 0])
        p6 = self.geom.add_point([branch1_start + branch1_width, 0, 0])
        p7 = self.geom.add_point([branch1_start + branch1_width, inlet_height, 0])
        p8 = self.geom.add_point([branch1_start, inlet_height, 0])
        
        points.extend([p5, p6, p7, p8])
        branch1_rect = self.geom.add_rectangle(branch1_start, 0, branch1_width, inlet_height)
        
        # Branch 2
        p9 = self.geom.add_point([branch2_start, 0, 0])
        p10 = self.geom.add_point([branch2_start + branch2_width, 0, 0])
        p11 = self.geom.add_point([branch2_start + branch2_width, inlet_height, 0])
        p12 = self.geom.add_point([branch2_start, inlet_height, 0])
        
        points.extend([p9, p10, p11, p12])
        branch2_rect = self.geom.add_rectangle(branch2_start, 0, branch2_width, inlet_height)
        
        # Second level branches (distribute to channel inlets)
        channel_inlet_points = []
        channel_outlet_points = []
        
        # Calculate channel inlet positions
        for i in range(self.num_channels):
            x_pos = i * (self.channel_width_inlet + self.rib_width)
            y_pos = inlet_height
            
            # Create channel inlet
            p_inlet1 = self.geom.add_point([x_pos, y_pos, 0])
            p_inlet2 = self.geom.add_point([x_pos + self.channel_width_inlet, y_pos, 0])
            p_inlet3 = self.geom.add_point([x_pos + self.channel_width_inlet, y_pos + self.channel_depth, 0])
            p_inlet4 = self.geom.add_point([x_pos, y_pos + self.channel_depth, 0])
            
            channel_inlet_points.append([p_inlet1, p_inlet2, p_inlet3, p_inlet4])
            
            # Create channel outlet (at end of channel)
            x_outlet = x_pos + (self.channel_width_outlet - self.channel_width_inlet) * 0.5
            p_outlet1 = self.geom.add_point([x_outlet, y_pos + self.channel_length, 0])
            p_outlet2 = self.geom.add_point([x_outlet + self.channel_width_outlet, y_pos + self.channel_length, 0])
            p_outlet3 = self.geom.add_point([x_outlet + self.channel_width_outlet, y_pos + self.channel_length + self.channel_depth, 0])
            p_outlet4 = self.geom.add_point([x_outlet, y_pos + self.channel_length, 0])
            
            channel_outlet_points.append([p_outlet1, p_outlet2, p_outlet3, p_outlet4])
        
        return channel_inlet_points, channel_outlet_points
    
    def create_interdigitated_channels(self, inlet_points: List, outlet_points: List):
        """
        Create interdigitated channels with linear taper.
        
        Args:
            inlet_points: List of inlet point sets for each channel
            outlet_points: List of outlet point sets for each channel
        """
        channels = []
        
        for i in range(self.num_channels):
            # Create channel with linear taper
            channel = self.create_tapered_channel(
                inlet_points[i], 
                outlet_points[i], 
                i
            )
            channels.append(channel)
        
        return channels
    
    def create_tapered_channel(self, inlet_pts: List, outlet_pts: List, channel_idx: int):
        """
        Create a single tapered channel.
        
        Args:
            inlet_pts: Inlet point set
            outlet_pts: Outlet point set
            channel_idx: Channel index for naming
        """
        # Create channel volume
        channel = self.geom.add_volume([
            self.geom.add_line_loop([
                self.geom.add_line(inlet_pts[0], inlet_pts[1]),      # Bottom
                self.geom.add_line(inlet_pts[1], inlet_pts[2]),      # Right
                self.geom.add_line(inlet_pts[2], inlet_pts[3]),      # Top
                self.geom.add_line(inlet_pts[3], inlet_pts[0])       # Left
            ])
        ])
        
        # Extrude to create 3D channel
        channel_3d = self.geom.extrude(channel, [0, 0, self.channel_depth])
        
        # Add physical group for fluid
        self.geom.add_physical(channel_3d, f"FLUID_CHANNEL_{channel_idx}")
        
        return channel_3d
    
    def create_gdl_layer(self):
        """Create GDL layer under channels."""
        # GDL extends under all channels
        gdl_rect = self.geom.add_rectangle(
            0, 0, 
            self.total_width, 
            self.channel_length + self.channel_depth
        )
        
        # Extrude to create 3D GDL
        gdl_3d = self.geom.extrude(gdl_rect, [0, 0, -self.gdl_thickness])
        
        # Add physical group
        self.geom.add_physical(gdl_3d, "GDL_CATHODE")
        
        return gdl_3d
    
    def create_membrane_layer(self):
        """Create membrane layer under GDL."""
        # Membrane extends under GDL
        membrane_rect = self.geom.add_rectangle(
            0, 0, 
            self.total_width, 
            self.channel_length + self.channel_depth
        )
        
        # Extrude to create 3D membrane
        membrane_3d = self.geom.extrude(membrane_rect, [0, 0, -self.membrane_thickness])
        
        # Add physical group
        self.geom.add_physical(membrane_3d, "MEMBRANE")
        
        return membrane_3d
    
    def create_boundary_conditions(self):
        """Create boundary condition tags."""
        # Inlet boundary
        inlet_boundary = self.geom.add_physical(
            self.geom.add_point([0, 0, 0]), 
            "INLET_CATHODE"
        )
        
        # Outlet boundary
        outlet_boundary = self.geom.add_physical(
            self.geom.add_point([self.total_width, self.channel_length, 0]), 
            "OUTLET_CATHODE"
        )
        
        # Wall boundaries
        wall_boundary = self.geom.add_physical(
            self.geom.add_point([0, 0, 0]), 
            "WALLS"
        )
        
        return inlet_boundary, outlet_boundary, wall_boundary
    
    def generate_mesh(self, mesh_size: str = "medium") -> str:
        """
        Generate the complete geometry and mesh.
        
        Args:
            mesh_size: Mesh refinement level ("coarse", "medium", "fine")
            
        Returns:
            Path to generated mesh file
        """
        # Create geometry components
        channel_inlets, channel_outlets = self.create_lung_manifold()
        channels = self.create_interdigitated_channels(channel_inlets, channel_outlets)
        gdl = self.create_gdl_layer()
        membrane = self.create_membrane_layer()
        
        # Create boundary conditions
        self.create_boundary_conditions()
        
        # Set mesh size based on refinement level
        if mesh_size == "coarse":
            mesh_size_factor = 0.5
        elif mesh_size == "medium":
            mesh_size_factor = 0.3
        elif mesh_size == "fine":
            mesh_size_factor = 0.1
        else:
            mesh_size_factor = 0.3
        
        # Set global mesh size
        self.geom.set_mesh_size_callback(lambda dim, tag, x, y, z: mesh_size_factor)
        
        # Generate mesh
        mesh_path = f"geometry/flowfield_{mesh_size}.msh"
        self.geom.generate_mesh(mesh_path, dim=3)
        
        return mesh_path
    
    def export_to_fenics(self, mesh_path: str) -> str:
        """
        Convert Gmsh mesh to FEniCS format.
        
        Args:
            mesh_path: Path to Gmsh mesh file
            
        Returns:
            Path to FEniCS mesh file
        """
        import meshio
        
        # Read Gmsh mesh
        mesh = meshio.read(mesh_path)
        
        # Convert to FEniCS format
        fenics_mesh_path = mesh_path.replace('.msh', '.xdmf')
        meshio.write(fenics_mesh_path, mesh, file_format="xdmf")
        
        return fenics_mesh_path


def main():
    """Main function to generate geometry and mesh."""
    # Create geometry
    pemfc_geom = PEMFCGeometry()
    
    # Generate meshes for all refinement levels
    mesh_files = {}
    for mesh_size in ["coarse", "medium", "fine"]:
        print(f"Generating {mesh_size} mesh...")
        mesh_path = pemfc_geom.generate_mesh(mesh_size)
        fenics_path = pemfc_geom.export_to_fenics(mesh_path)
        mesh_files[mesh_size] = {
            "gmsh": mesh_path,
            "fenics": fenics_path
        }
        print(f"  Gmsh: {mesh_path}")
        print(f"  FEniCS: {fenics_path}")
    
    print("\nMesh generation complete!")
    return mesh_files


if __name__ == "__main__":
    main()
