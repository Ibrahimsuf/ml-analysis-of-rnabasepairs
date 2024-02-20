import pandas as pd
from preprocess_data import Preprocessor
import numpy as np
import os
from tqdm import tqdm
import json

class GetInteratomicDistances:
    def __init__(self, pdbs_dir, base_pair_distances_file) -> None:
        self.pdbs_dir = pdbs_dir
        self.base_pair_distances_file = base_pair_distances_file
        self.base_pair_distances_df = pd.read_csv(self.base_pair_distances_file)
        # self.distances = {"AA": [],"AU": [], "AC": [], "AG": [], "CC": [], "CG": [], "CU": [], "GG": [], "GU": [], "UU": []}

        G_atoms = ['P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'O3\'', 'C2\'', 'O2\'', 'C1\'', 'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4']
        A_atoms = ['P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', "O4'", 'C3\'', 'O3\'', 'C2\'', 'O2\'', 'C1\'', 'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4']
        C_atoms = ['P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'O3\'', 'C2\'', 'O2\'', 'C1\'', 'N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6']
        U_atoms = ['P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'O3\'', 'C2\'', 'O2\'', 'C1\'', 'N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6']

        self.atoms_in_bases = {"G": G_atoms, "A": A_atoms, "C": C_atoms, "U": U_atoms}
        self.incomplete_bases = []
    def get_all_interatomic_distances(self, distances_path, errors_path):
        pdb_list = self.base_pair_distances_df["pdb_id"].unique().tolist()

        for pdb_id in pdb_list[:1]:
            print("Finding interatomic distances for {}".format(pdb_id))
            distances = self.get_interatomic_distances_for_pdb(pdb_id)

            for basepair_type, data in distances.items():
                pd.DataFrame(data).to_csv(os.path.join(distances_path, f"{basepair_type}.csv"), mode="w", index=False, header=True)

        for pdb_id in tqdm(pdb_list[1:]):
            print("Finding interatomic distances for {}".format(pdb_id))
            distances = self.get_interatomic_distances_for_pdb(pdb_id)

            for basepair_type, data in distances.items():
                pd.DataFrame(data).to_csv(os.path.join(distances_path, f"{basepair_type}.csv"), mode="a", index=False, header=False)




        pd.DataFrame(self.incomplete_bases).to_csv(errors_path)
        
    def get_interatomic_distances_for_pdb(self, pdb_id):
        distances = {"AA": [],"AU": [], "AC": [], "AG": [], "CC": [], "CG": [], "CU": [], "GG": [], "GU": [], "UU": []}
        pdb_possible_pairs = self.base_pair_distances_df[self.base_pair_distances_df["pdb_id"].str.replace(".pdb", "").str.replace(".cif", "") == pdb_id.replace(".pdb", "").replace(".cif", "")]
        model = Preprocessor.read_pdb(pdb_id, self.pdbs_dir)
        for _, basepair_row in pdb_possible_pairs.iterrows():
            record = {}
            record["pdb_id"] = pdb_id
            record["nt1"] = basepair_row["nt1"]
            record["nt2"] = basepair_row["nt2"]
            record["BasePair"] = basepair_row["BasePair"]

            residue1 = model[str(basepair_row["chain1"])][eval(basepair_row["residue1"])]
            residue2 = model[str(basepair_row["chain2"])][eval(basepair_row["residue2"])]

            res1name = Preprocessor.replaceHetatms(residue1)
            res2name = Preprocessor.replaceHetatms(residue2)

            chain1 = basepair_row["chain1"]
            chain2 = basepair_row["chain2"]

            if not res1name in ["A", "C", "G", "U"] or not res2name in ["A", "C", "G", "U"]:
                raise(NameError) 
            

            # We switch the base pairs if the later in the alphabet base pair comes first so the distances are consistent, but keep the name so it matches with dssr
            if (res1name + res2name) in distances:
                basepair_type = res1name + res2name
            elif (res2name + res1name) in distances:
                basepair_type = res2name + res1name
                residue1, residue2 = residue2, residue1
                chain1, chain2 = chain2, chain1
                res1name, res2name = res2name, res1name
            else:
                print(basepair_row["chain1"])
                print(basepair_row["residue1"])
                print(basepair_row["chain2"])
                print(basepair_row["residue2"])
                print(res2name + res1name)
                raise(NameError)

            #for now we are skipping A1 and D1 because they are missing some of the atoms in the residue eventually we will need to find a way to deal with this
            #and missing values in general
            
            requied_atoms_in_base1 = set(self.atoms_in_bases[res1name])
            atom_names_in_residue1 = set([atom.get_name() for atom in residue1])

            requied_atoms_in_base2 = set(self.atoms_in_bases[res2name])
            atom_names_in_residue2 = set([atom.get_name() for atom in residue2])


            if not requied_atoms_in_base1.issubset(atom_names_in_residue1):
                self.incomplete_bases.append((pdb_id, basepair_row["chain1"], eval(basepair_row["residue1"]), tuple(requied_atoms_in_base1 - atom_names_in_residue1)))
                continue
            if not requied_atoms_in_base2.issubset(atom_names_in_residue2):
                self.incomplete_bases.append((pdb_id, basepair_row["chain2"], eval(basepair_row["residue2"]), tuple(requied_atoms_in_base2 - atom_names_in_residue2)))
                continue


            for atom1 in requied_atoms_in_base1:
                for atom2 in requied_atoms_in_base2:
                    # we convert to float 64 so it is json serializeable
                    record[f"{atom1}-{atom2}"] = np.float64(residue1[atom1] - residue2[atom2])
            distances[basepair_type].append(record)

        return distances
    



def main():
    # getinteratomicdistances = GetInteratomicDistances("/Users/ibrahims/Documents/Programming/undergrad_reasearch/rna/rna_parsing/pdbs_0_3", "0_3/residue_pair_distances.csv")
    # getinteratomicdistances.get_all_interatomic_distances("0_3/interatomic_distances.json", "0_3/errors.csv")

    # getinteratomicdistances = GetInteratomicDistances("/Users/ibrahims/Documents/Programming/undergrad_reasearch/rna/rna_annotations_parser/pdbs/3_5/cif_files/", "3_5/residue_pair_distances.csv")
    # getinteratomicdistances.get_all_interatomic_distances("3_5", "3_5/errors.csv")
    
    GetInteratomicDistances("test_pdbs", "test_pdbs/residue_pair_distances.csv").get_all_interatomic_distances("test_pdbs", "test_pdbs/errors.csv")


if __name__ == "__main__":
    main()