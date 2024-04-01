import csv
import pandas as pd
import os
import gzip
class AddStacksandPairstoDistances:
    def __init__(self, stacks_dir, pairs_dir, base_pair_distances_file):
        self.stacks_dir = stacks_dir
        self.parirs_dir = pairs_dir
        self.base_pair_distances_file = base_pair_distances_file
        self.base_pair_distances = pd.read_csv(base_pair_distances_file)
        self.pdb_id = None
        self.stacks_file = None
        self.base_stacks = None
        self.base_pairs_file = None
    def addstacksandpairs(self):
      with gzip.open(self.base_pair_distances_file,'rt') as csvinput:
        with gzip.open(f'{self.base_pair_distances_file[:-4]}_withlabels.csv.', 'wt') as csvoutput:
            writer = csv.writer(csvoutput, lineterminator='\n')
            reader = csv.reader(csvinput)
            row = next(reader)
            row.append("Class")
            writer.writerow(row)

            for row in reader:
                pdb_id, chain1, chain2, residue1, residue2 = row[25:30]
                if pdb_id != self.pdb_id:
                    #if we are on to a new pdb load that into memory and create a list of the base stacks
                    self.pdb_id = pdb_id
                    print(self.pdb_id)
                    self.basepairs_file = pd.read_csv(f"{self.parirs_dir}/{pdb_id.upper().replace('.PDB', '').replace('.CIF', '')}.csv", index_col=0)
                    self.basestacks_file = pd.read_csv(f"{self.stacks_dir}/{pdb_id.upper().replace('.PDB', '').replace('.CIF', '')}.csv")
                    self.basepairs_file["nt1"] = self.basepairs_file["nt1"].str.replace("/", "") if "nt1" in self.basepairs_file.columns else -1
                    self.basepairs_file["nt2"] = self.basepairs_file["nt2"].str.replace("/", "") if "nt2" in self.basepairs_file.columns else -1
                    self.basepairs_file.set_index(["nt1", "nt2"], inplace=True)
                    self.basestacks_file = self.basestacks_file.set_index(["nt1", "nt2"])
                
                nt1 = f"{chain1}{eval(residue1)[1]}"
                nt2 = f"{chain2}{eval(residue2)[1]}"
                row.append(self.getclass(nt1, nt2))
                writer.writerow(row)

    def getclass(self, nt1, nt2):
      if (nt1, nt2) in self.basepairs_file.index or (nt2, nt1) in self.basepairs_file.index:
        _class = 1
      else:
        if (nt1, nt2) in self.basestacks_file.index or (nt2, nt1) in self.basestacks_file.index:
          _class = 2
        else:
          _class = 0
      return _class


def main():
    resolution = "3_5"
    basepairtypes = ["UU"]
    for basepairtype in basepairtypes:
        stacks_dir = f"dssrstacks_{resolution}"
        pairs_dir = f"dssrpairs{resolution}"
        base_pair_distances_file = f"{resolution}/{basepairtype}_beaddistances.csv.gz"
        AddStacksandPairstoDistances(stacks_dir, pairs_dir, base_pair_distances_file).addstacksandpairs()

        # os.remove(base_pair_distances_file)

    # for basepairtype in basepairtypes:
    #     stacks_dir = "dssrstacks_3_5"
    #     base_pair_distances_file = f"3_5/{basepairtype}.csv"
    #     AddBaseStackstoDistances(stacks_dir, base_pair_distances_file, index=False).addstackstofile()

    os.remove(base_pair_distances_file)


if __name__ == "__main__":
    main()