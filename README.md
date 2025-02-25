# RNALig (Scoring function to predict RNA-ligand binding affinity)
RNALig is an AI/ML-based scoring function designed to predict RNA-ligand binding affinity. By leveraging advanced machine learning algorithms, RNALig analyzes the interactions between RNA molecules and ligands, providing a quantitative prediction of their binding strength. This tool helps streamline the discovery of RNA-targeting drugs and aids in understanding RNA-ligand binding mechanisms. It integrates various features, including structural data and chemical properties, to create accurate models. The function's high efficiency and predictive power make it a valuable resource for researchers in drug design, molecular biology, and bioinformatics, facilitating the identification of potential therapeutic agents targeting RNA.

**Steps to use**
**Affinity conversion:** AI/ML algorithms work on same unit. So to train the model all the experimental binding affinities (kd) has been converted into kcal/mol to train the model.
**Features extraction:** The next step to extract the essential features. So the features has been divided into three parts: RNAspecific features, Ligand specific features and Complex specific features. the Features_extraction.py script will be used for this.
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
After features extraction use 










