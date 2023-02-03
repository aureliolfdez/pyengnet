from setuptools import setup, find_packages
import os
from os import path

this_directory = path.abspath(path.dirname(__file__))

# read the contents of README.rst
def readme():
    with open(path.join(this_directory, 'README.rst'), encoding='utf-8') as f:
        return f.read()

# Setting up
setup(
    name="pyengnet",
    version='0.0.3',
    author="Aurelio Lopez-Fernandez",
    author_email="alopfer1@upo.es",
    description='pyEnGNet: optimized reconstruction of gene co-expression networks using multi-GPU',    
    long_description=readme(),
    long_description_content_type="text/x-rst",
    packages=['pyengnet'],
    install_requires=['numpy', 'tqdm', 'multiprocess', 'pandas', 'matplotlib', 'networkx', 'scipy'],
    keywords=['python', 'multiGPU', 'bioinformatics', 'gene networks'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering'
    ],    
    python_requires='>=3.9'
)