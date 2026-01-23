from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
this_directory = Path(__file__).parent
long_description = (
    (this_directory / "README.md").read_text()
    if (this_directory / "README.md").exists()
    else ""
)

setup(
    name="orthoSynAssign",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Ortholog Synteny Assignment Tool for analyzing OrthoFinder results with genome annotations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/orthoSynAssign",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        # Add any dependencies here
        # Example: "pandas>=1.0.0",
    ],
    entry_points={
        "console_scripts": [
            "orthoSynAssign=orthoSynAssign:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
