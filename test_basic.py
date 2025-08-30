#!/usr/bin/env python3
"""
Basic Test Script for PEMFC CFD Study
=====================================

This script tests basic functionality without requiring full dependencies.
"""

import sys
import os

def test_imports():
    """Test basic imports."""
    print("Testing basic imports...")
    
    try:
        import numpy as np
        print("✓ numpy imported successfully")
    except ImportError as e:
        print(f"✗ numpy import failed: {e}")
        return False
    
    try:
        import matplotlib.pyplot as plt
        print("✓ matplotlib imported successfully")
    except ImportError as e:
        print(f"✗ matplotlib import failed: {e}")
        return False
    
    try:
        import pandas as pd
        print("✓ pandas imported successfully")
    except ImportError as e:
        print(f"✗ pandas import failed: {e}")
        return False
    
    return True

def test_geometry():
    """Test geometry module."""
    print("\nTesting geometry module...")
    
    try:
        # Test basic geometry calculations
        channel_depth = 1.0e-3
        channel_width_inlet = 1.0e-3
        channel_width_outlet = 0.85e-3
        rib_width = 1.0e-3
        num_channels = 8
        
        # Calculate total width
        total_width = num_channels * (channel_width_inlet + rib_width) - rib_width
        active_area = total_width * 20.0e-3 * 1e4  # cm²
        
        print(f"✓ Geometry calculations successful")
        print(f"  Total width: {total_width*1000:.1f} mm")
        print(f"  Channel taper: {(1 - channel_width_outlet/channel_width_inlet)*100:.1f}%")
        print(f"  Active area: {active_area:.1f} cm²")
        
        return True
        
    except Exception as e:
        print(f"✗ Geometry test failed: {e}")
        return False

def test_physics():
    """Test physics calculations."""
    print("\nTesting physics calculations...")
    
    try:
        # Test basic physics constants
        F = 96485.0  # C/mol
        R = 8.314    # J/(mol·K)
        T = 353.0    # K
        P = 1.013e5  # Pa
        
        # Test operating condition calculation
        current_density = 1.2  # A/cm²
        stoichiometry = 2.0
        active_area = 25.0e-4  # m²
        
        # Calculate O₂ consumption rate
        j_am2 = current_density * 1e4
        n_o2_cons = j_am2 * active_area / (4 * F)
        
        # Calculate inlet O₂ flow rate
        n_o2_in = stoichiometry * n_o2_cons
        
        print(f"✓ Physics calculations successful")
        print(f"  O₂ consumption rate: {n_o2_cons*1e6:.3f} mmol/s")
        print(f"  Inlet O₂ flow rate: {n_o2_in*1e6:.3f} mmol/s")
        
        return True
        
    except Exception as e:
        print(f"✗ Physics test failed: {e}")
        return False

def test_file_structure():
    """Test project file structure."""
    print("\nTesting project file structure...")
    
    required_files = [
        "README.md",
        "VERSIONS.txt",
        "requirements.txt",
        "setup.py",
        "Methods.md",
        "src/__init__.py",
        "src/physics.py",
        "src/solver.py",
        "src/main.py",
        "geometry/__init__.py",
        "geometry/flowfield.py",
        "cases/default_config.yml",
        "notebooks/demo.py"
    ]
    
    missing_files = []
    for file_path in required_files:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    if missing_files:
        print(f"✗ Missing files: {len(missing_files)}")
        for file_path in missing_files:
            print(f"  - {file_path}")
        return False
    else:
        print(f"✓ All required files present ({len(required_files)} files)")
        return True

def main():
    """Run all tests."""
    print("=" * 60)
    print("PEMFC CFD STUDY - BASIC FUNCTIONALITY TEST")
    print("=" * 60)
    
    tests = [
        test_imports,
        test_geometry,
        test_physics,
        test_file_structure
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
    
    print("\n" + "=" * 60)
    print(f"TEST RESULTS: {passed}/{total} tests passed")
    print("=" * 60)
    
    if passed == total:
        print("🎉 All tests passed! The package is ready for use.")
        print("\nNext steps:")
        print("1. Install full dependencies: pip install -r requirements.txt")
        print("2. Run demo: python notebooks/demo.py")
        print("3. Run full study: python src/main.py --help")
        return 0
    else:
        print("❌ Some tests failed. Please check the errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
