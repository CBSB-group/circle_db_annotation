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

## Instructions

### 1. Download Project Files
Download and extract the `variant_annotation_pipeline` directory. Ensure the following files are present:

* `variant_annotation.py` – Main pipeline script
* `mutate_model.py` – Adapted MODELLER script for *in silico* mutations
* `misc.py` – Helper functions for the pipeline
* `chimera_script.sh` – Shell script for Chimera 5 Å environments generation
* `P06280.pdb` – Example PDB structure file
* `P06280.csv` – Example variants input file

---

### 2. Configure System PATH
Add `pdb2pqr`, `apbs`, and `chimera` to your system's `PATH` environment variable:

```bash
echo 'export PATH="/path/to/pdb2pqr/bin:$PATH"' >> ~/.bashrc
echo 'export PATH="/path/to/apbs/bin:$PATH"' >> ~/.bashrc
echo 'export PATH="/path/to/chimera/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```
### 3. Test

Execute the main script using the provided example input files (`P06280.csv` and `P06280.pdb`):

```bash
cd variant_annotation
python variant_annotation.py -i P06280.csv -p P06280.pdb
```
### 📂 Generated Output

Execution creates the **`variant_annotation_output.csv`** dataset.

> [!NOTE]
> All parameters listed below are reported in pairs: for both the **Wild-Type** residue and the corresponding **Mutated** residue.

* **Secondary structure** (DSSP mapping)
* **Solvent accessibility** (SASA calculation)
* **Charge assignment** (PDB2PQR generated)
* **Electrostatic potential** (APBS calculation)
