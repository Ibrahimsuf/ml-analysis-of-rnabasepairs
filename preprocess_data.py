from Bio.PDB import PDBParser
import os
import numpy as np
import pandas as pd
import re
from tqdm import tqdm

class Preprocessor:
    def __init__(self, pdbs_dir, annotations_folder, cutoff) -> None:
        self.pdbs_dir = pdbs_dir
        self.annotations_folder = annotations_folder
        self.residue_pair_distances_df = pd.DataFrame(columns = ["pdb_id", "chain1", "chain2", "residue1", "residue2", "distance", "nt1", "nt2"])
        self.cutoff = cutoff

    def preproces_pdbs(self):
        pdb_list = os.listdir(self.pdbs_dir)
        for pdb_id in tqdm(pdb_list):
            self._preprocess_pdb(pdb_id)

    def _preprocess_pdb(self, pdb_id):
        residue_pair_distances = self.get_base_pairs_within_cutoff(pdb_id)
        dssr_labels = self.get_dssr_annotations(pdb_id)
        pairs_df = residue_pair_distances.merge(dssr_labels, on = ["nt1", "nt2"], how="outer")

        if pairs_df["distance"].isnull().any():
            raise Exception("PDB {} has missing distances".format(pdb_id))
        
        self.residue_pair_distances_df = pd.concat([self.residue_pair_distances_df, pairs_df])


    def _read_pdb(self, pdb_id):
        # pdb_id = pdb_id.upper() + ".pdb"
        # print(pdb_id)
        parser = PDBParser()
        path = os.path.join(self.pdbs_dir, pdb_id)
        structure =  parser.get_structure(pdb_id, path)[0]
        return structure
    
    def get_base_pairs_within_cutoff(self, pdb_id):
        model = self._read_pdb(pdb_id)
        residue_pair_distances = {"pdb_id": [], "chain1": [], "chain2": [], "residue1": [], "residue2": [], "distance": []}
       
        chain_names = []
        for chain in model:
            chain_names.append(chain.id)

        for chain1 in tqdm(model):

            for chain2 in tqdm(model):
                if chain_names.index(chain2.id) > chain_names.index(chain2.id):
                    continue
                for residue1 in tqdm(chain1):
                    res1name = Preprocessor.replaceHetatms(residue1)
                    if not res1name in ["A", "U", "G", "C"]:
                        continue
                    for residue2 in chain2:
                        res2name = Preprocessor.replaceHetatms(residue2)
                        if not res2name in ["A", "U", "G", "C"]:
                            continue
                        if chain1.id == chain2.id and residue1.id[1] > residue2.id[1]:
                            continue
                        if residue1 != residue2:
                            center1 = residue1.center_of_mass(geometric=True)
                            center2 = residue2.center_of_mass(geometric=True)
                            distance = np.linalg.norm(center1 - center2)

                            if distance > self.cutoff:
                                continue
                            residue_pair_distances["pdb_id"].append(pdb_id.replace(".pdb", ""))
                            residue_pair_distances["chain1"].append(chain1.id)
                            residue_pair_distances["chain2"].append(chain2.id)
                            residue_pair_distances["residue1"].append(residue1.id)
                            residue_pair_distances["residue2"].append(residue2.id)

                            residue_pair_distances["distance"].append(distance)
        
        # print(pdb_id)
        # print(residue_pair_distances)
        residue_pair_distances = pd.DataFrame(residue_pair_distances)
        # print(residue_pair_distances.apply(lambda x: x["chain1"] + str(x["residue1"][1]), axis = 1))
        residue_pair_distances["nt1"] = residue_pair_distances.apply(lambda x: x["chain1"] + str(x["residue1"][1]), axis = 1)
        residue_pair_distances["nt2"] = residue_pair_distances.apply(lambda x: x["chain2"] + str(x["residue2"][1]), axis = 1)

        return residue_pair_distances
    
    def get_dssr_annotations(self, pdb_id):
        annotations_file = os.path.join(self.annotations_folder, pdb_id.lower().replace(".pdb", ""), "dssr_annotations.csv")
        dssr_labels = pd.read_csv(annotations_file)
        if len(dssr_labels) == 0:
            return pd.DataFrame(columns = ["nt1", "nt2", "BasePair"])
        
        dssr_labels["nt1"] = dssr_labels["nt1"].apply(lambda x: x.replace("/", ""))
        dssr_labels["nt2"] = dssr_labels["nt2"].apply(lambda x: x.replace("/", ""))
        bases_to_filter_out = ['P', "p", 'a', 'u', "c", 'g', "T", "T"] # Lowercase bases are modified bases we are ignoring these for now


        pattern = '|'.join(map(re.escape, bases_to_filter_out))
        dssr_labels = dssr_labels[~dssr_labels["bp"].str.contains(pattern, case=True)]

        dssr_labels = dssr_labels[["nt1", 'nt2']]
        dssr_labels["BasePair"] = True

        dssr_labels = dssr_labels.apply(Preprocessor.swap_residues, axis=1)
        
        return dssr_labels
    
    @staticmethod
    def swap_residues(dsssr_label_row):
        residue1_number = ''.join(re.findall(r'\d', dsssr_label_row["nt1"]))
        residue2_number = ''.join(re.findall(r'\d', dsssr_label_row["nt2"]))

        residue1_letter = ''.join(re.findall(r'[a-zA-Z]', dsssr_label_row["nt1"]))
        residue2_letter = ''.join(re.findall(r'[a-zA-Z]', dsssr_label_row["nt2"]))
        
        if int(residue1_number) > int(residue2_number) and residue1_number > residue2_number and residue1_letter == residue2_letter:
            dsssr_label_row["nt1"], dsssr_label_row["nt2"] = dsssr_label_row["nt2"], dsssr_label_row["nt1"]
        
        return dsssr_label_row
    
    @staticmethod
    def replaceHetatms(residue):
        # we are removing the basepairs that have modified bases in them instead of chainging them into the normal bases


        # mappings = {"A2M": "A", "IU": "U", "CCC": "C", "GTP": "G", "C2E": "G", "ADE": "A", "1MA": "A", "OMU": "U"}
        mappings = {"ADE": "A", "DC": "C", "DA": "A", "DU": "U", "DG": "G"}
        resname = residue.resname
        for key, val in mappings.items():
            resname = resname.replace(key, val)

        return resname


def main():
    preprocessor = Preprocessor("/Users/ibrahims/Documents/Programming/undergrad_reasearch/rna/rna_parsing/pdbs_0_3", "/Users/ibrahims/Documents/Programming/undergrad_reasearch/rna/rna_parsing", 15)
    preprocessor.preproces_pdbs()


if __name__ == "__main__":
    main()