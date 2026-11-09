import subprocess
import os
from Bio.PDB import PDBParser, DSSP


def delete_file(percorso):
    if os.path.exists(percorso):
        if os.path.isfile(percorso):
            os.remove(percorso)
            print(f"File deleted: {percorso}")
        else:
            print(f"The path is not a file: {percorso}")
    else:
        print(f"File not found: {percorso}")


def get_residue_one_letter(pdb_file: str, chain_id: str, residue_number: int):
    """
    Returns the one-letter code of the residue found
    at a given position in the PDB.
    """

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]

    chain = model[chain_id]

    for residue in chain:

        # Ignore HETATM and water
        if residue.id[0] != " ":
            continue

        if residue.id[1] == residue_number:

            resname = residue.get_resname().upper()

            if resname not in AA3_TO_AA1:
                raise ValueError(
                    f"Unrecognized residue: {resname}"
                )

            return AA3_TO_AA1[resname]

    raise ValueError(
        f"Residue {residue_number} not found in chain {chain_id}."
    )



def run_pdb2pqr_and_apbs(
    pdb_path: str,
    forcefield: str,
    output_pqr_path: str,
    chain_id: str,
    ph: float,
    apbs_input_path: str,
    apbs_output_file: str,
    residue_position: int
):

    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB not found: {pdb_path}")

    if not chain_id or len(chain_id) != 1:
        raise ValueError("chain_id must be a single character (e.g. 'A')")

    # Run pdb2pqr
    cmd_pdb2pqr = [
        "pdb2pqr",
        f"--ff={forcefield}",
        "--chain", chain_id,
        f"--with-ph={ph}",
        "--apbs-input", apbs_input_path,
        pdb_path,
        output_pqr_path
    ]

    try:
        subprocess.run(cmd_pdb2pqr, check=True)

        log_path = os.path.splitext(output_pqr_path)[0] + ".log"
        if os.path.exists(log_path):
            os.remove(log_path)

        print("pdb2pqr completed.")

    except subprocess.CalledProcessError as e:
        raise RuntimeError("pdb2pqr failed") from e

    # Run APBS
    try:
        subprocess.run([
            "apbs",
            apbs_input_path,
            "--output-format=flat",
            f"--output-file={apbs_output_file}"
        ], check=True)

        print("APBS completed.")

    except subprocess.CalledProcessError as e:
        raise RuntimeError("APBS failed") from e

    # Read PQR: residue charges
    residue_charges = {}

    with open(output_pqr_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                parts = line.split()
                res_number = int(parts[4])
                charge = float(parts[8])

                residue_charges.setdefault(res_number, 0.0)
                residue_charges[res_number] += charge

    if residue_position not in residue_charges:
        raise ValueError(
            f"Residue {residue_position} not found in the PQR file"
        )

    net_charge = round(residue_charges[residue_position], 4)

    # Read energy from the .pot file
    pot = None
    inside_calc = False
    current_id = None

    with open(apbs_output_file, "r") as f:
        for line in f:
            stripped = line.strip()

            if stripped.startswith("calc"):
                inside_calc = True
                current_id = None

            elif inside_calc and stripped.startswith("id"):
                current_id = int(stripped.split()[1])

            elif inside_calc and stripped.startswith("totEnergy") and current_id == 2:
                pot = round(float(stripped.split()[1]), 3)

            elif stripped.startswith("end"):
                inside_calc = False

    if pot is None:
        raise ValueError("totEnergy for id 2 not found in the .pot file")

    pot = round(pot, 4)

    print(f"Residue {residue_position} net charge: {net_charge}")
    print(f"TotEnergy id 2: {pot} kJ/mol")

    # Remove temporary files
    files_to_remove = [
        output_pqr_path,
        output_pqr_path + ".dx",
        apbs_output_file,
        apbs_input_path,
        "io.mc"
    ]

    for fpath in files_to_remove:
        if os.path.exists(fpath):
            os.remove(fpath)

    print("Temporary files removed.")

    return net_charge, pot
    

def get_secondary_structure(pdb_file: str, chain_id: str, residue_number: int, return_type: str = "ss"):
    """
    Returns DSSP information for a specific residue in a PDB file.

    Args:
        pdb_file (str): path to the PDB file
        chain_id (str): chain ID (e.g. 'A')
        residue_number (int): residue number in the PDB
        return_type (str): 'ss' for secondary structure, 'acc' for solvent accessibility

    Returns:
        tuple: (residue_position, value)
            residue_position = (chain_id, residue_number)
            value = DSSP code or relative accessibility
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]  # first model

    # Compute DSSP
    dssp = DSSP(model, pdb_file)

    # DSSP uses keys like: (chain_id, (' ', resseq, icode))
    key = (chain_id, (' ', residue_number, ' '))

    if key not in dssp:
        raise ValueError(
            f"Residue {residue_number} in chain {chain_id} not found in DSSP."
        )

    ss = dssp[key][2]      # secondary structure code
    acc = dssp[key][3]     # relative accessibility

    if return_type == "ss":
        value = ss
    elif return_type == "acc":
        value = acc
    else:
        raise ValueError("return_type must be 'ss' or 'acc'")

    return value



def run_chimera_zone_script(pdb_file, chain, position, distance, output_file):
    try:
        subprocess.run(
            [
                "./chimera_script.sh",
                str(pdb_file),
                str(chain),
                str(position),
                str(distance),
                str(output_file)
            ],
            check=True
        )

        print(
            f"Residues within {distance} angstroms from {position}.{chain} "
            f"written to {output_file}."
        )

    except subprocess.CalledProcessError as e:
        print(f"Error while running Chimera: {str(e)}")


def mutate_residue(
    pdb_code: str,
    position: int,
    chain: str,
    aa_one_letter: str,
    log_file: str = "log.log"
):

    # One-letter -> three-letter dictionary (local to the function)
    aa_one_to_three = {
        "A": "ALA",
        "R": "ARG",
        "N": "ASN",
        "D": "ASP",
        "C": "CYS",
        "E": "GLU",
        "Q": "GLN",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "L": "LEU",
        "K": "LYS",
        "M": "MET",
        "F": "PHE",
        "P": "PRO",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
        "Y": "TYR",
        "V": "VAL",
    }

    aa_one_letter = aa_one_letter.upper()

    if aa_one_letter not in aa_one_to_three:
        raise ValueError(f"Invalid amino acid: {aa_one_letter}")

    aa_three_letter = aa_one_to_three[aa_one_letter]

    # Command with shell redirection
    command = (
        f"python mutate_model.py "
        f"{pdb_code} {position} {aa_three_letter} {chain} "
        f"> {log_file}"
    )

    # Run the command
    subprocess.run(command, shell=True, check=True)

    # Delete the log file
    if os.path.exists(log_file):
        os.remove(log_file)
        
        
        
AA3_TO_AA1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


