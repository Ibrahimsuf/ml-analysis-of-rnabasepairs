import csv
import pandas as pd
import os
class AddBaseStackstoDistances:
    def __init__(self, stacks_dir, base_pair_distances_file, index = True):
        self.stacks_dir = stacks_dir
        self.base_pair_distances_file = base_pair_distances_file
        self.base_pair_distances = pd.read_csv(base_pair_distances_file)
        self.pdb_id = None
        self.stacks_file = None
        self.base_stacks = None
        self.index = index

    def addstackstofile(self):
        with open(self.base_pair_distances_file,'r') as csvinput:
            with open(f'{self.base_pair_distances_file[:-4]}_stacks{self.base_pair_distances_file[-4:]}', 'w') as csvoutput:
                writer = csv.writer(csvoutput, lineterminator='\n')
                reader = csv.reader(csvinput)
                row = next(reader)
                row.append('BaseStack')
                writer.writerow(row)

                for row in reader:
                    if self.index:
                        nt1, nt2 = row[2], row[3]
                    else:
                        nt1, nt2 = row[1], row[1]
                    
                    if self.index:
                        pdb_id = row[1]
                    else:
                        pdb_id = row[0]
                    if pdb_id != self.pdb_id:
                        #if we are on to a new pdb load that into memory and create a list of the base stacks
                        self.pdb_id = pdb_id
                        print(self.pdb_id)
                        self.stacks_file = pd.read_csv(f"{self.stacks_dir}/{self.pdb_id.replace('.pdb', '').replace('.cif', '')}.csv")
                        self.stacks_file["residues_sorted"] = self.stacks_file.apply(lambda row: tuple(sorted([str(row["nt1"]), str(row["nt2"])])), axis = 1)
                        self.base_stacks = set(self.stacks_file["residues_sorted"].values)
                    
                    row.append(self.findifresiduesarestack(nt1, nt2))
                    writer.writerow(row)

    def findifresiduesarestack(self, nt1, nt2):
        residue_pair = tuple(sorted([nt1, nt2]))
        if residue_pair in self.base_stacks:
            return True
        else:
            return False


def main():
    basepairtypes = ["AA", "AC", "AG", "AU", "CC", "CG", "CU", "GG", "GU", "UU"]

    for basepairtype in basepairtypes:
        stacks_dir = "dssrstacks_0_3"
        base_pair_distances_file = f"0_3/{basepairtype}.csv"
        AddBaseStackstoDistances(stacks_dir, base_pair_distances_file).addstackstofile()

        os.remove(base_pair_distances_file)

    for basepairtype in basepairtypes:
        stacks_dir = "dssrstacks_3_5"
        base_pair_distances_file = f"3_5/{basepairtype}.csv"
        AddBaseStackstoDistances(stacks_dir, base_pair_distances_file, index=False).addstackstofile()

        os.remove(base_pair_distances_file)


if __name__ == "__main__":
    main()