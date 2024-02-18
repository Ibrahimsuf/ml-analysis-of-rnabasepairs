import warnings
import pandas as pd
import os
from preprocess_data import Preprocessor
warnings.filterwarnings('ignore')

def getcoords(pdbs_dir, base_pair_distances_df):
    pdbs = os.listdir(pdbs_dir)
    coords = []

    for pdb_id in pdbs:
        pdb_base_possible_pair = base_pair_distances_df[base_pair_distances_df["pdb_id"].str.replace(".pdb", "").str.replace(".cif", "") == pdb_id.replace(".pdb", "").replace(".cif", "")]
        model = Preprocessor.read_pdb(pdb_id)
        print(pdb_id)
        for _, basepair_row in pdb_base_possible_pair.iterrows():
            record = {}
            record["pdb_id"] = pdb_id
            record["nt1"] = basepair_row["nt1"]
            record["nt2"] = basepair_row["nt2"]
            record["BasePair"] = basepair_row["BasePair"]
            # We add the stacks latter because they are not in the distances file
            residue1 = model[basepair_row["chain1"]][eval(basepair_row["residue1"])]
            residue2 = model[basepair_row["chain2"]][eval(basepair_row["residue2"])]

            res1name = Preprocessor.replaceHetatms(residue1)
            res2name = Preprocessor.replaceHetatms(residue2)
            if not res1name in ["A", "C", "G", "U"] or not res2name in ["A", "C", "G", "U"]:
                raise(NameError) 
        
            for atom in residue1:
                record[f"r1-{atom.get_name()}"] = tuple(atom.get_coord())
            for atom in residue2:
                record[f"r2-{atom.get_name()}"] = tuple(atom.get_coord())
            coords.append(record)
    coords_df = pd.DataFrame(coords)
    return coords_df


def main():
    pdbs_dir = "/Users/ibrahims/Documents/Programming/undergrad_reasearch/rna/rna_parsing/pdbs_0_3"
    base_pair_distances_df = pd.read_csv("0_3/residue_pair_distances.csv")
    coords_df = getcoords(pdbs_dir, base_pair_distances_df)
    coords_df.to_csv("0_3/coords.csv")


    pdbs_dir = "/Users/ibrahims/Documents/Programming/undergrad_reasearch/rna/rna_annotations_parser/pdbs/3_5/cif_files"
    base_pair_distances_df = pd.read_csv("3_5/residue_pair_distances.csv")
    coords_df = getcoords(pdbs_dir, base_pair_distances_df)
    coords_df.to_csv("3_5/coords.csv")