import warnings
import pandas as pd
import os
import csv
from preprocess_data import Preprocessor
warnings.filterwarnings('ignore')


def get_column_names(pdbs_dir,base_pair_distances_df):
    column_names = {"pdb_id", "nt1", "nt2", "BasePair"}
    pdbs = os.listdir(pdbs_dir)
    for pdb_id in pdbs:
        pdb_base_possible_pair = base_pair_distances_df[base_pair_distances_df["pdb_id"].str.replace(".pdb", "").str.replace(".cif", "") == pdb_id.replace(".pdb", "").replace(".cif", "")]
        model = Preprocessor.read_pdb(pdb_id, pdbs_dir=pdbs_dir)
        print(pdb_id)
        for _, basepair_row in pdb_base_possible_pair.iterrows():
            # We add the stacks latter because they are not in the distances file
            residue1 = model[str(basepair_row["chain1"])][eval(basepair_row["residue1"])]
            residue2 = model[str(basepair_row["chain2"])][eval(basepair_row["residue2"])]

            res1name = Preprocessor.replaceHetatms(residue1)
            res2name = Preprocessor.replaceHetatms(residue2)
            if not res1name in ["A", "C", "G", "U"] or not res2name in ["A", "C", "G", "U"]:
                raise(NameError) 
        
            for atom in residue1:
                column_names.add(f"r1-{atom.get_name()}")
            for atom in residue2:
                column_names.add(f"r2-{atom.get_name()}")
    return column_names
def getcoords(pdbs_dir, base_pair_distances_df, outfile):
    column_names = get_column_names(pdbs_dir, base_pair_distances_df)
    print(column_names)
    pdbs = os.listdir(pdbs_dir)
    with open(outfile, "w") as f:
        writer = csv.DictWriter(f, fieldnames=column_names)
        writer.writeheader()
        for pdb_id in pdbs:
            pdb_base_possible_pair = base_pair_distances_df[base_pair_distances_df["pdb_id"].str.replace(".pdb", "").str.replace(".cif", "") == pdb_id.replace(".pdb", "").replace(".cif", "")]
            model = Preprocessor.read_pdb(pdb_id, pdbs_dir=pdbs_dir)
            for _, basepair_row in pdb_base_possible_pair.iterrows():
                record = {}
                record["pdb_id"] = pdb_id
                record["nt1"] = basepair_row["nt1"]
                record["nt2"] = basepair_row["nt2"]
                record["BasePair"] = basepair_row["BasePair"]
                # We add the stacks latter because they are not in the distances file
                residue1 = model[str(basepair_row["chain1"])][eval(basepair_row["residue1"])]
                residue2 = model[str(basepair_row["chain2"])][eval(basepair_row["residue2"])]

                res1name = Preprocessor.replaceHetatms(residue1)
                res2name = Preprocessor.replaceHetatms(residue2)
                if not res1name in ["A", "C", "G", "U"] or not res2name in ["A", "C", "G", "U"]:
                    raise(NameError) 
            
                for atom in residue1:
                    record[f"r1-{atom.get_name()}"] = tuple(atom.get_coord())
                for atom in residue2:
                    record[f"r2-{atom.get_name()}"] = tuple(atom.get_coord())

                writer.writerow(record)


def main():
    # pdbs_dir = "/Users/ibrahims/Documents/Programming/undergrad_reasearch/rna/rna_parsing/pdbs_0_3"
    # base_pair_distances_df = pd.read_csv("0_3/residue_pair_distances.csv")
    # coords_df = getcoords(pdbs_dir, base_pair_distances_df)
    # coords_df.to_csv("0_3/coords.csv")
    pdbs_dir = "/Users/ibrahims/Documents/Programming/undergrad_reasearch/rna/rna_annotations_parser/pdbs/3_5/cif_files"
    base_pair_distances_df = pd.read_csv("3_5/residue_pair_distances.csv")
    getcoords(pdbs_dir, base_pair_distances_df, "3_5/coords.csv")

if __name__ == "__main__":
    main()