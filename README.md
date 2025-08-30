# PEMFC Bio-Inspired CFD Study

## Overview
This repository contains a comprehensive CFD study of a PEM fuel cell with a **tapered interdigitated cathode flow field** fed by a **two-level lung-inspired inlet manifold**. The study uses FEniCSx (dolfinx) with PETSc for high-performance computing and Gmsh for mesh generation.

## Key Features
- **Geometry**: Tapered interdigitated channels (10-20% area reduction downstream)
- **Manifold**: Two-level lung-inspired branching (area-preserving ratios)
- **Physics**: Multi-physics coupling (flow, heat, species, electrochemistry)
- **Validation**: Mesh independence study with three refinement levels
- **Output**: Publication-ready figures and comprehensive data analysis

## Quick Start (Google Colab)
```bash
# Clone repository
!git clone https://github.com/yourusername/pemfc-bio-inspired-cfd.git
cd pemfc-bio-inspired-cfd

# Install dependencies
!pip install fenics-dolfinx petsc4py mpi4py gmsh pygmsh meshio numpy scipy pandas matplotlib pyvista

# Run complete study
!python src/main.py --mesh fine --current 1.2 --stoich high
```

## Local Setup
```bash
# Clone and setup
git clone https://github.com/yourusername/pemfc-bio-inspired-cfd.git
cd pemfc-bio-inspired-cfd

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run study
python src/main.py --mesh fine --current 1.2 --stoich high
```

## Project Structure
```
pemfc-fenics/
├── geometry/              # pygmsh sources + exported .msh
├── meshes/                # coarse/medium/fine meshes
├── src/                   # solvers & tools
├── cases/                 # YAML configs
├── outputs/
│   ├── figures/          # PNG + PDF
│   ├── data/             # CSV tables
│   └── logs/             # solver logs
├── notebooks/             # demo notebooks
├── Methods.md             # detailed methodology
├── README.md              # this file
└── VERSIONS.txt          # environment specs
```

## Operating Conditions
- **Current Densities**: 0.3, 0.8, 1.2 A/cm²
- **Stoichiometries**: Low (1.5) and High (2.5)
- **Temperature**: 353 K (80°C)
- **Pressure**: Atmospheric (1 atm)

## Outputs
- **Figures**: Velocity, temperature, species, and pressure fields
- **Data**: KPI tables, mesh independence reports, convergence histories
- **Analysis**: Performance comparison with conventional designs

## Citation
If you use this work in your research, please cite:
```
@software{pemfc_cfd_2024,
  title={PEMFC Bio-Inspired CFD Study},
  author={Your Name},
  year={2024},
  url={https://github.com/yourusername/pemfc-bio-inspired-cfd}
}
```

## License
MIT License - see LICENSE file for details.
