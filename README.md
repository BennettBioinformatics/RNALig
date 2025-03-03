# RNALig (Scoring function to predict RNA-ligand binding affinity)
RNALig is an AI/ML-based scoring function designed to predict RNA-ligand binding affinity. By leveraging advanced machine learning algorithms, RNALig analyzes the interactions between RNA molecules and ligands, providing a quantitative prediction of their binding strength. This tool helps streamline the discovery of RNA-targeting drugs and aids in understanding RNA-ligand binding mechanisms. It integrates various features, including structural data and chemical properties, to create accurate models. The function's high efficiency and predictive power make it a valuable resource for researchers in drug design, molecular biology, and bioinformatics, facilitating the identification of potential therapeutic agents targeting RNA.

**Steps to use**
**Features extraction:** The first is step to extract the essential features. So the features has been divided into three parts: RNAspecific features, Ligand specific features and Complex specific features. the Features_extraction.py script will be used for this.
import os
import csv
import MDAnalysis as mda
from Bio.PDB import PDBParser
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Lipinski, Crippen, rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
import subprocess
import numpy as np
from scipy.spatial.distance import cdist
These essential packages need to be installed before use the script.

**Binding affinity prediction**
After features extraction use Random Forest Refressor script for binding affinity prediction. 
Before run the script put the features file in working folder and install all the required libraries. The script will save the binding affinities in the form of .csv file. The script is able to predict binding affinity for one RNA-ligand complex or many number of complexes. So during running of script it will popup to ask for specific pdb structure or number of pdbs. So if user have one pdb structure or several structures should be save in working folder only. So script directly took the structure from working folder and predict affinity.
Both features extraction and Random Forest Regressor will work on same working folder.

**Affinity Conversion**
Affinity conversion script has been used to convert the experimental affinities in one unit kcal/mol to train the model. If user wants to convert the affinity then this script will be helpful.

**Install packages**
1- MDAnalysis
#pip install MDAnalysis
@varify the installation (python or python3)
#python3 -c "import MDAnalysis; print(MDAnalysis.__version__)"
@Install other dependencies
#pip install biopython pubchempy rdkit scipy numpy
@Install ViennaRNA 
sudo apt-get install build-essential libgsl-dev libgslcblas0 zlib1g-dev
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_5_x/ViennaRNA-2.5.0.tar.gz
tar -xvzf ViennaRNA-2.5.0.tar.gz
cd ViennaRNA-2.5.0
./configure --prefix=/usr/local
make
sudo make install
After installing these packages the Features_extraction.py script will work properly.
2- 











