import os
import argparse
import pandas as pd

# Helper functions live in misc.py: they wrap the external tools
# (Modeller, DSSP, Chimera, PDB2PQR + APBS) behind simple Python calls.
from misc import (
    mutate_residue,
    get_secondary_structure,
    get_residue_one_letter,
    run_chimera_zone_script,
    run_pdb2pqr_and_apbs
)

# Force field and pH used by PDB2PQR when assigning charges / protonation
# states. They are global so that WT and mutant are always processed with
# identical settings, which keeps the two sets of values comparable.
forcefield="PARSE"
ph=7
# Used when the user does not provide -O, or provides only a directory.
DEFAULT_OUTPUT_NAME = "variant_annotation_output.csv"

def parse_args():
    # Command-line interface: -i and -p are mandatory, -o and -O are optional.
    parser = argparse.ArgumentParser(
        description="Variant annotation on a PDB structure."
    )

    # Variants table; required columns are validated later in main().
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="CSV file containing the variants."
    )

    # Structure on which the mutations are modelled and annotated.
    parser.add_argument(
        "-p",
        "--pdb",
        required=True,
        help="PDB file to use."
    )

    # Lets the user map Uniprot numbering onto the residue numbering of the
    # PDB file, which often differs (e.g. when the construct is truncated).
    parser.add_argument(
        "-o",
        "--offset",
        type=int,
        default=0,
        help="Offset to add to the Uniprot numbering (default: 0)."
    )

    # Note the case: lowercase -o is the offset, uppercase -O is the output.
    parser.add_argument(
        "-O",
        "--output",
        default=os.path.join(os.getcwd(), DEFAULT_OUTPUT_NAME),
        help=(
            "Output CSV file name and/or path. If a directory is given, "
            f"'{DEFAULT_OUTPUT_NAME}' is created inside it "
            "(default: current directory/%(default)s)."
        )
    )

    return parser.parse_args()


def main():

    args = parse_args()

    # Resolve the output path
    # (expanduser handles "~", abspath makes it independent of the cwd)
    output_path = os.path.abspath(os.path.expanduser(args.output))

    # If the user passed a directory, use the default file name inside it
    # The endswith check is needed because abspath strips the trailing slash,
    # so a not-yet-existing directory such as "results/" would otherwise be
    # mistaken for a file name.
    if os.path.isdir(output_path) or args.output.endswith(("/", os.sep)):
        output_path = os.path.join(output_path, DEFAULT_OUTPUT_NAME)

    # Create the output directory if it does not exist
    # (done now so the pipeline does not fail only at the final to_csv call,
    # after all the expensive computations have already been run)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Read the CSV
    df = pd.read_csv(args.input)

    # Check required columns
    # Failing early with a clear message is better than a KeyError in the
    # middle of the loop.
    required_columns = [
        "Uniprot",
        "PDB",
        "Chain",
        "OrigAA",
        "Position",
        "SubAA",
    ]

    missing = [c for c in required_columns if c not in df.columns]

    if missing:
        raise ValueError(
            f"Missing required columns: {', '.join(missing)}"
        )

    # Output columns
    # Created up front (empty) so every row ends up with the same schema,
    # whether it succeeds or is later flagged as "Missing".
    # Orig_* = wild-type residue, Sub_* = mutated residue.
    df["MutatedFile"] = ""
    df["Orig_Sec"] = ""
    df["Sub_Sec"] = ""
    df["Orig_Acc"] = ""
    df["Sub_Acc"] = ""
    df["Orig_Charge"] = ""
    df["Sub_Charge"] = ""
    df["Orig_Pot"] = ""
    df["Sub_Pot"] = ""

    # Echo the run configuration for traceability in logs.
    print(f"CSV    : {args.input}")
    print(f"PDB    : {args.pdb}")
    print(f"Offset : {args.offset}")
    print(f"Output : {output_path}")
    print()

    # Loop over the variants
    # Each variant is processed independently: an error on one row must not
    # stop the others (see the try/except below).
    for idx, row in df.iterrows():

        uniprot = row["Uniprot"]
        pdb_name = row["PDB"]
        chain = row["Chain"]

        orig_aa = row["OrigAA"]
        position = int(row["Position"])
        sub_aa = row["SubAA"]

        # Apply the offset
        # From here on, "position" is the residue number in the PDB file,
        # not the Uniprot one, and it is used for every downstream call.
        position = position + args.offset

        # idx is 0-based, so +1 gives a human-friendly progress counter.
        print(
            f"[{idx + 1}/{len(df)}] "
            f"{orig_aa}{position}{sub_aa}"
        )

        # Name of the mutant structure; it must match the file name that
        # mutate_residue produces, since it is read back by name below.
        mutated_file = f"{uniprot}_{position}{sub_aa}.pdb"

        try:
            # Initialised before anything can fail, so the finally block can
            # safely check them even if the error occurs early.
            wt_5angfile = None
            mut_5angfile = None            

            # Check that the WT residue matches
            # Guards against a wrong offset or a wrong PDB/chain: modelling
            # a mutation on the wrong residue would silently give garbage.
            pdb_residue = get_residue_one_letter(
                pdb_file=pdb_name,
                chain_id=chain,
                residue_number=position
            )

            if pdb_residue != orig_aa:
                raise ValueError(
                    f"Expected residue {orig_aa}, found {pdb_residue}"
                )

            # Generate the mutant
            # In silico mutagenesis with Modeller; writes mutated_file to disk.
            mutate_residue(
                pdb_code=pdb_name,
                position=position,
                chain=chain,
                aa_one_letter=sub_aa
            )

            # Recorded only after the mutation succeeded.
            df.at[idx, "MutatedFile"] = mutated_file

            # ==========================
            # Wild type
            # ==========================

            # DSSP secondary structure of the WT residue.
            df.at[idx, "Orig_Sec"] = get_secondary_structure(
                pdb_file=pdb_name,
                chain_id=chain,
                residue_number=position,
                return_type="ss"
            )

            # DSSP solvent accessibility of the WT residue, rounded to
            # 3 decimals to keep the CSV compact.
            df.at[idx, "Orig_Acc"] = round(
                get_secondary_structure(
                    pdb_file=pdb_name,
                    chain_id=chain,
                    residue_number=position,
                    return_type="acc"
                ),
                3
            )

            # ==========================
            # Mutant
            # ==========================

            # Same DSSP annotations, now on the modelled mutant structure,
            # so that WT and mutant values can be compared directly.
            df.at[idx, "Sub_Sec"] = get_secondary_structure(
                pdb_file=mutated_file,
                chain_id=chain,
                residue_number=position,
                return_type="ss"
            )

            df.at[idx, "Sub_Acc"] = round(
                get_secondary_structure(
                    pdb_file=mutated_file,
                    chain_id=chain,
                    residue_number=position,
                    return_type="acc"
                ),
                3
            )
            
            # Local environment, WT: Chimera extracts the atoms within 5 Å of
            # the residue of interest. Electrostatics are then computed on this
            # small environment only, which keeps APBS fast and focuses the
            # result on the immediate surroundings of the variant.
            wt_5angfile = pdb_name.replace(".pdb", "_5ang.pdb")

            run_chimera_zone_script(
                pdb_file=pdb_name,
                chain=chain,
                position=position,
                distance=5,
                output_file=wt_5angfile
            )
            
            # Local environment, mutant: same 5 Å selection on the modelled
            # structure, to compare like with like.
            mut_5angfile = mutated_file.replace(".pdb", "_5ang.pdb")

            run_chimera_zone_script(
                pdb_file=mutated_file,
                chain=chain,
                position=position,
                distance=5,
                output_file=mut_5angfile
            )
            
            # Intermediate files for the WT electrostatics, all derived from
            # the 5 Å file name: PQR (charges/radii from PDB2PQR), APBS input
            # file, and the potential output file.
            wt_pqr_path = wt_5angfile.replace(".pdb", ".pqr")
            wt_apbs_input = wt_5angfile.replace(".pdb", ".in")
            wt_apbs_output = wt_5angfile.replace(".pdb", ".pot")
            


            # PDB2PQR assigns charges/protonation (forcefield and pH set at the
            # top of the file) and generates the APBS input; APBS then computes
            # the potential. Returns the net charge and the electrostatic
            # potential, the latter evaluated at the variant position.
            wt_net_charge, wt_pot = run_pdb2pqr_and_apbs(
                pdb_path=wt_5angfile,
                forcefield=forcefield,
                output_pqr_path=wt_pqr_path,
                chain_id=chain,
                ph=ph,
                apbs_input_path=wt_apbs_input,
                apbs_output_file=wt_apbs_output,
                residue_position=position
            )

            df.at[idx, "Orig_Charge"] = wt_net_charge
            df.at[idx, "Orig_Pot"] = wt_pot
            
            # Same intermediate file names for the mutant environment.
            mut_pqr_path = mut_5angfile.replace(".pdb", ".pqr")
            mut_apbs_input = mut_5angfile.replace(".pdb", ".in")
            mut_apbs_output = mut_5angfile.replace(".pdb", ".pot")
            

            # Same charge + electrostatics calculation on the mutant, with
            # identical force field and pH as the WT run.
            mut_net_charge, mut_pot = run_pdb2pqr_and_apbs(
                pdb_path=mut_5angfile,
                forcefield=forcefield,
                output_pqr_path=mut_pqr_path,
                chain_id=chain,
                ph=ph,
                apbs_input_path=mut_apbs_input,
                apbs_output_file=mut_apbs_output,
                residue_position=position
            )

            df.at[idx, "Sub_Charge"] = mut_net_charge
            df.at[idx, "Sub_Pot"] = mut_pot
            

        except Exception as e:

            # Any failure (residue mismatch, Modeller/Chimera/APBS error, ...)
            # is reported and the row is marked "Missing" in every output
            # column, so partially filled values are never mixed with failures.
            print(f"Error at row {idx + 1}: {e}")

            df.at[idx, "MutatedFile"] = "Missing"
            df.at[idx, "Orig_Sec"] = "Missing"
            df.at[idx, "Sub_Sec"] = "Missing"
            df.at[idx, "Orig_Acc"] = "Missing"
            df.at[idx, "Sub_Acc"] = "Missing"
            df.at[idx, "Orig_Charge"] = "Missing"
            df.at[idx, "Sub_Charge"] = "Missing"
            df.at[idx, "Orig_Pot"] = "Missing"
            df.at[idx, "Sub_Pot"] = "Missing"

        finally:

            # Cleanup runs on both success and failure, so that temporary
            # structures do not pile up in the working directory.
            # Delete the mutated file if it exists
            if os.path.exists(mutated_file):
                os.remove(mutated_file)

            # The None checks cover the case where the error happened before
            # these paths were assigned.
            if wt_5angfile and os.path.exists(wt_5angfile):
                os.remove(wt_5angfile)

            if mut_5angfile and os.path.exists(mut_5angfile):
                os.remove(mut_5angfile)

    # Save the output CSV
    # Written once at the end, after all variants have been processed.
    df.to_csv(output_path, index=False)

    print(f"\nPipeline completed. Output saved to: {output_path}")


if __name__ == "__main__":
    main()
