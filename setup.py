from setuptools import setup, find_packages
import os

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

# Setting up
setup(
    name="pyengnet",
    version='0.0.1.4',
    author="Aurelio Lopez-Fernandez",
    author_email="alopfer1@upo.es",
    description='pyEnGNet: optimized reconstruction of gene co-expression networks using multi-GPU',
    long_description_content_type="text/markdown",
    long_description='pyEnGNet: optimized reconstruction of gene co-expression networks using multi-GPU',
    package_dir={'': 'src'},
    packages=find_packages("src"),
    install_requires=['numpy', 'tqdm', 'multiprocess', 'pandas', 'matplotlib', 'networkx', 'scipy'],
    keywords=['python', 'multiGPU', 'bioinformatics', 'vgene networks'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10'    
)