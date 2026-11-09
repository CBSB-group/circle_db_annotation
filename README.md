# circle_db_annotation

A Python-based pipeline for the characterization of missense protein variants.

# Overview

This pipeline automates the annotation of protein missense variants. Given a PDB structure and a csv file of mutations, it performs in silico mutagenesis, generates local structural environments, and calculates secondary structure, solvent accessibility, charges, and electrostatic potentials.

## Dependencies and Requirements

- **Python** `v3.11.8` – Core programming language.
- **Biopython** – Handling and parsing of PDB files.
- **Modeller** `v10.6` – *In silico* mutations.
- **Chimera** `v1.19` *(build 42556)* – Generation of 5 Å environments.
- **DSSP / mkdssp** `v3.0.0` – Secondary structure and accessibility assignment.
- **PDB2PQR** `v3.6.2` – Generation of charges and APBS input files.
- **APBS** `v3.4.1` – Electrostatic potential calculation.

