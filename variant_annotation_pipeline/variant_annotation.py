#!/usr/bin/env python3

import os
import argparse
import pandas as pd

from misc import (
    mutate_residue,
    get_secondary_structure,
    get_residue_one_letter,
    run_chimera_zone_script,
    run_pdb2pqr_and_apbs
)

forcefield="PARSE"
ph=7
def parse_args():
    parser = argparse.ArgumentParser(
        description="Variant annotation on a PDB structure."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="CSV file containing the variants."
    )

    parser.add_argument(
        "-p",
        "--pdb",
        required=True,
        help="PDB file to use."
    )

    parser.add_argument(
        "-o",
        "--offset",
        type=int,
        default=0,
        help="Offset to add to the Uniprot numbering (default: 0)."
    )

    return parser.parse_args()


def main():

    args = parse_args()

    # Read the CSV
    df = pd.read_csv(args.input)

    # Check required columns
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
    df["MutatedFile"] = ""
    df["Orig_Sec"] = ""
    df["Sub_Sec"] = ""
    df["Orig_Acc"] = ""
    df["Sub_Acc"] = ""
    df["Orig_Charge"] = ""
    df["Sub_Charge"] = ""
    df["Orig_Pot"] = ""
    df["Sub_Pot"] = ""

    print(f"CSV    : {args.input}")
    print(f"PDB    : {args.pdb}")
    print(f"Offset : {args.offset}")
    print()

    # Loop over the variants
    for idx, row in df.iterrows():

        uniprot = row["Uniprot"]
        pdb_name = row["PDB"]
        chain = row["Chain"]

        orig_aa = row["OrigAA"]
        position = int(row["Position"])
        sub_aa = row["SubAA"]

        # Apply the offset
        position = position + args.offset

        print(
            f"[{idx + 1}/{len(df)}] "
            f"{orig_aa}{position}{sub_aa}"
        )

        mutated_file = f"{uniprot}_{position}{sub_aa}.pdb"

        try:
            wt_5angfile = None
            mut_5angfile = None            

            # Check that the WT residue matches
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
            mutate_residue(
                pdb_code=pdb_name,
                position=position,
                chain=chain,
                aa_one_letter=sub_aa
            )

            df.at[idx, "MutatedFile"] = mutated_file

            # ==========================
            # Wild type
            # ==========================

            df.at[idx, "Orig_Sec"] = get_secondary_structure(
                pdb_file=pdb_name,
                chain_id=chain,
                residue_number=position,
                return_type="ss"
            )

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
            
            wt_5angfile = pdb_name.replace(".pdb", "_5ang.pdb")

            run_chimera_zone_script(
                pdb_file=pdb_name,
                chain=chain,
                position=position,
                distance=5,
                output_file=wt_5angfile
            )
            
            mut_5angfile = mutated_file.replace(".pdb", "_5ang.pdb")

            run_chimera_zone_script(
                pdb_file=mutated_file,
                chain=chain,
                position=position,
                distance=5,
                output_file=mut_5angfile
            )
            
            wt_pqr_path = wt_5angfile.replace(".pdb", ".pqr")
            wt_apbs_input = wt_5angfile.replace(".pdb", ".in")
            wt_apbs_output = wt_5angfile.replace(".pdb", ".pot")
            


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
            
            mut_pqr_path = mut_5angfile.replace(".pdb", ".pqr")
            mut_apbs_input = mut_5angfile.replace(".pdb", ".in")
            mut_apbs_output = mut_5angfile.replace(".pdb", ".pot")
            

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

            # Delete the mutated file if it exists
            if os.path.exists(mutated_file):
                os.remove(mutated_file)

            if wt_5angfile and os.path.exists(wt_5angfile):
                os.remove(wt_5angfile)

            if mut_5angfile and os.path.exists(mut_5angfile):
                os.remove(mut_5angfile)

    # Save the output CSV
    df.to_csv("variant_annotation_output.csv", index=False)

    print("\nPipeline completed.")


if __name__ == "__main__":
    main()
