#!/usr/bin/env python3
"""
Setup script for PEMFC CFD Study Package
"""

from setuptools import setup, find_packages
import os

# Read requirements
with open('requirements.txt', 'r') as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]

# Read README
with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="pemfc-cfd-study",
    version="1.0.0",
    description="PEMFC Bio-Inspired CFD Study using FEniCSx",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="PEMFC Research Team",
    author_email="research@pemfc-study.org",
    url="https://github.com/yourusername/pemfc-bio-inspired-cfd",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.9",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
            "mypy>=0.800",
        ],
        "docs": [
            "sphinx>=4.0",
            "sphinx-rtd-theme>=1.0",
            "nbsphinx>=0.8",
        ],
    },
    entry_points={
        "console_scripts": [
            "pemfc-cfd=src.main:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
