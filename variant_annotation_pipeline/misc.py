import subprocess
import os
# Biopython: PDBParser reads the structure, DSSP is a wrapper around the
# external mkdssp executable (which must be installed and reachable).
from Bio.PDB import PDBParser, DSSP


# Small utility to remove a file with a friendly message. It is not used by
# variant_annotation.py.
def delete_file(file):
    if os.path.exists(file):
        # Refuse to delete anything that is not a regular file (e.g. a
        # directory), so the helper cannot remove more than intended.
        if os.path.isfile(file):
            os.remove(file)
            print(f"File deleted: {file}")
        else:
            print(f"The path is not a file: {file}")
    else:
        print(f"File not found: {file}")


# Used by the pipeline as a sanity check: the residue found in the structure
# at the given position must match the wild-type amino acid declared in the
# variants CSV, otherwise the offset, chain or PDB file is probably wrong.
def get_residue_one_letter(pdb_file: str, chain_id: str, residue_number: int):
    """
    Returns the one-letter code of the residue found
    at a given position in the PDB.
    """

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    # Only the first model is considered (relevant for NMR ensembles).
    model = structure[0]

    chain = model[chain_id]

    for residue in chain:

        # Ignore HETATM and water
        # In Biopython the first element of residue.id is the hetero flag:
        # " " for standard residues, "H_xxx" for hetero groups, "W" for water.
        if residue.id[0] != " ":
            continue

        # residue.id[1] is the residue number as written in the PDB file.
        if residue.id[1] == residue_number:

            resname = residue.get_resname().upper()

            # AA3_TO_AA1 is defined at the bottom of this module; it is
            # resolved at call time, so the order of definition does not matter.
            # Only the 20 standard amino acids are mapped: modified residues
            # (e.g. MSE) end up here and raise an error.
            if resname not in AA3_TO_AA1:
                raise ValueError(
                    f"Unrecognized residue: {resname}"
                )

            return AA3_TO_AA1[resname]

    raise ValueError(
        f"Residue {residue_number} not found in chain {chain_id}."
    )



# Assigns charges with PDB2PQR, runs APBS on the resulting input file and
# returns two values: the net charge of the residue at residue_position
# (summed over its atoms, from the PQR file) and the APBS total energy
# (kJ/mol) read from the output file. Intermediate files are cleaned up at
# the end of a successful run.
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

    # Fail early with a clear message: this is also what surfaces when the
    # Chimera step did not produce the 5 A environment file (see
    # run_chimera_zone_script, which does not raise on failure).
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB not found: {pdb_path}")

    if not chain_id or len(chain_id) != 1:
        raise ValueError("chain_id must be a single character (e.g. 'A')")

    # Run pdb2pqr
    # Besides the PQR file (charges and radii), --apbs-input makes pdb2pqr
    # write the input file that APBS needs in the next step.
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

        # pdb2pqr leaves a .log file next to the PQR; it is not needed.
        log_path = os.path.splitext(output_pqr_path)[0] + ".log"
        if os.path.exists(log_path):
            os.remove(log_path)

        print("pdb2pqr completed.")

    # check=True raises CalledProcessError on a non-zero exit code; it is
    # re-raised as RuntimeError so the caller's generic handler can report it.
    except subprocess.CalledProcessError as e:
        raise RuntimeError("pdb2pqr failed") from e

    # Run APBS
    # "flat" output format makes APBS write a plain-text results file, which
    # is parsed further down to extract the energy.
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
    # Sum of the atomic charges for every residue, keyed by residue number.
    residue_charges = {}

    with open(output_pqr_path, "r") as f:
        for line in f:
            # Only ATOM records carry charges (HETATM lines are skipped).
            if line.startswith("ATOM"):
                parts = line.split()
                # Positions in the whitespace-split record: 4 = residue
                # number, 8 = atomic charge. This depends on the exact column
                # layout of the PQR file written by pdb2pqr.
                res_number = int(parts[4])
                charge = float(parts[8])

                residue_charges.setdefault(res_number, 0.0)
                residue_charges[res_number] += charge

    # The residue of interest must be present in the 5 A environment file.
    if residue_position not in residue_charges:
        raise ValueError(
            f"Residue {residue_position} not found in the PQR file"
        )

    # Rounding removes floating-point noise from the sum of atomic charges.
    net_charge = round(residue_charges[residue_position], 4)

    # Read energy from the .pot file
    # Small state machine over the APBS output: the file is organised in
    # "calc ... end" blocks, each with an "id" and a "totEnergy" line.
    pot = None
    inside_calc = False
    current_id = None

    with open(apbs_output_file, "r") as f:
        for line in f:
            stripped = line.strip()

            # A new calculation block starts: reset the id.
            if stripped.startswith("calc"):
                inside_calc = True
                current_id = None

            elif inside_calc and stripped.startswith("id"):
                current_id = int(stripped.split()[1])

            # Only the totEnergy of the block with id 2 is kept; which
            # calculation this is depends on the APBS input generated by
            # pdb2pqr.
            elif inside_calc and stripped.startswith("totEnergy") and current_id == 2:
                pot = round(float(stripped.split()[1]), 3)

            elif stripped.startswith("end"):
                inside_calc = False

    # If the block was never found the value would silently stay None, so
    # fail explicitly.
    if pot is None:
        raise ValueError("totEnergy for id 2 not found in the .pot file")

    # Already rounded to 3 decimals above; this second rounding changes nothing.
    pot = round(pot, 4)

    print(f"Residue {residue_position} net charge: {net_charge}")
    print(f"TotEnergy id 2: {pot} kJ/mol")

    # Remove temporary files
    # Only reached when everything succeeded: after an error the intermediate
    # files are left in place, which can help with debugging. "io.mc" is the
    # log APBS writes in the working directory.
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
    

# Runs DSSP on the structure and returns either the secondary-structure code
# or the relative solvent accessibility of one residue. Note that, unlike
# what the docstring says, only the value is returned, not a tuple.
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
    # The PDB file is passed too because mkdssp reads it directly from disk.
    dssp = DSSP(model, pdb_file)

    # DSSP uses keys like: (chain_id, (' ', resseq, icode))
    # The key is built with an empty insertion code (' '), so residues with
    # an insertion code cannot be found this way.
    key = (chain_id, (' ', residue_number, ' '))

    # Residues that DSSP skips (e.g. incomplete backbone) are not in the
    # result and would otherwise cause a KeyError.
    if key not in dssp:
        raise ValueError(
            f"Residue {residue_number} in chain {chain_id} not found in DSSP."
        )

    # Fields of the DSSP record: index 2 = secondary structure, index 3 =
    # relative accessible surface area.
    ss = dssp[key][2]      # secondary structure code
    acc = dssp[key][3]     # relative accessibility

    if return_type == "ss":
        value = ss
    elif return_type == "acc":
        value = acc
    else:
        raise ValueError("return_type must be 'ss' or 'acc'")

    return value



# Calls the Chimera wrapper script to write the atoms within "distance"
# angstroms of the given residue to output_file. The script is looked up in
# the current directory ("./"), so the pipeline must be run from the folder
# that contains chimera_script.sh.
def run_chimera_zone_script(pdb_file, chain, position, distance, output_file):
    try:
        subprocess.run(
            [
                "./chimera_script.sh",
                # Arguments are cast to str because subprocess only accepts
                # strings, while position and distance are integers.
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

    # Unlike the other helpers, a failure here is only printed and not
    # re-raised: the pipeline then fails at the next step, when the expected
    # output file turns out to be missing.
    except subprocess.CalledProcessError as e:
        print(f"Error while running Chimera: {str(e)}")


# Creates the mutant structure by running the Modeller-based script
# mutate_model.py. The script receives the structure (pdb_code, which the
# pipeline fills with the PDB file name), the position, the three-letter code
# of the new amino acid and the chain.
def mutate_residue(
    pdb_code: str,
    position: int,
    chain: str,
    aa_one_letter: str,
    log_file: str = "log.log"
):

    # One-letter -> three-letter dictionary (local to the function)
    # Modeller expects three-letter residue names, while the CSV uses
    # one-letter codes.
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

    # Accept lowercase input too.
    aa_one_letter = aa_one_letter.upper()

    if aa_one_letter not in aa_one_to_three:
        raise ValueError(f"Invalid amino acid: {aa_one_letter}")

    aa_three_letter = aa_one_to_three[aa_one_letter]

    # Command with shell redirection
    # Modeller is very verbose: its output is redirected to a log file to
    # keep the console readable.
    command = (
        f"python3 mutate_model.py "
        f"{pdb_code} {position} {aa_three_letter} {chain} "
        f"> {log_file}"
    )

    # Run the command
    # shell=True is needed for the ">" redirection; check=True raises if the
    # script fails, in which case the log file is kept (the deletion below is
    # not reached), so it can be inspected.
    subprocess.run(command, shell=True, check=True)

    # Delete the log file
    if os.path.exists(log_file):
        os.remove(log_file)
        
        
        
# Lookup table used by get_residue_one_letter. It is defined at the end of the
# module but is only needed when that function is called.
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
