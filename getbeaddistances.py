import pandas as pd
import numpy as np
import csv
import argparse

def getdistances(beadchain1, beadchain2, beadsinresidueforrbead1, beadsinresidueforrbead2):
  distances = {}
  for bead1 in beadsinresidueforrbead1:
    for bead2 in beadsinresidueforrbead2:
      distances[f"{bead1}-{bead2}"] = np.linalg.norm(np.array(eval(beadchain1[f"{bead1}_coords"])) - np.array(eval(beadchain2[f"{bead2}_coords"])))

  return distances

def main(residue_pair_distances_path, beads_path, output_path):
  residue_pair_distances = pd.read_csv(residue_pair_distances_path)
  beads = pd.read_csv(beads_path, index_col=0)
  beads = beads.set_index(["pdb_id", "chain", "residue"])
  beadsinresidue = {"A": ["P", "S", "R1", "A1", "A2"], "U": ["P", "S", "Y1", "Y2", "U1"], "G": ["P", "S", "R1", "G1", "G2"], "C": ["P", "S", "Y1", "Y2", "C1"]}

  basepairfiles = {}
  for basepairtype in ["AA", "AC", "AG", "AU", "CC", "CG", "CU", "GG", "GU", "UU"]:
    basepairfiles[basepairtype] = open(f"{output_path}/{basepairtype}_beaddistances.csv", "w")


  for _, row in residue_pair_distances.iterrows():
    beadchain1 = beads.loc[row["pdb_id"], row["chain1"], row["residue1"]] if (row["pdb_id"], row["chain1"], row["residue1"]) in beads.index else None
    beadchain2 = beads.loc[row["pdb_id"], row["chain2"], row["residue2"]] if (row["pdb_id"], row["chain2"], row["residue2"]) in beads.index else None
    if beadchain1 is None or beadchain2 is None:
      continue
    
    basepair = "".join(sorted(beadchain1['resname']+beadchain2['resname']))
    beadsinresidueforrbead1 = beadsinresidue[beadchain1["resname"]]
    beadsinresidueforrbead2 = beadsinresidue[beadchain2["resname"]]

    if beadchain1["resname"] > beadchain2["resname"]:
      beadchain1, beadchain2 = beadchain2, beadchain1
      beadsinresidueforrbead1, beadsinresidueforrbead2 = beadsinresidueforrbead2, beadsinresidueforrbead1


    distances = getdistances(beadchain1, beadchain2, beadsinresidueforrbead1, beadsinresidueforrbead2)
    distances["pdb_id"] = row["pdb_id"]
    distances["chain1"] = row["chain1"]
    distances["chain2"] = row["chain2"]
    distances["residue1"] = row["residue1"]
    distances["residue2"] = row["residue2"]

    file = basepairfiles[basepair]
    writer = csv.DictWriter(file, distances.keys())
    if file.tell() == 0:
      writer.writeheader()
    writer.writerow(distances)

  

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("residue_pair_distances_path")
  parser.add_argument("beads_path")
  parser.add_argument("output_path")
  args = parser.parse_args()
  main(args.residue_pair_distances_path, args.beads_path, args.output_path)

