import csv
import pandas as pd
import os



def addclassificationtodata():
  for basepairtype in ["AA", "AC", "AG", "AU", "CC", "CG", "CU", "GG", "GU", "UU"]:
    addclassificationtodata_for_basepairtype(basepairtype)
    os.remove(f"0_3/{basepairtype}_stacks.csv")

def addclassificationtodata_for_basepairtype(basepairtype):
  with open(f"0_3/{basepairtype}_stacks.csv", "r") as input_file, open(f"0_3/{basepairtype}_stacks_classifcation.csv", "w") as output_file:
    reader = csv.reader(input_file)
    writer = csv.writer(output_file)
    nt1_index = 2
    nt2_index = 3
    pdb_id_index = 1
    basepair_index = 4
    current_pdb_id = None
    current_dssr_annotations = None
    row1 = next(reader)
    writer.writerow(row1 + ["Classification"])
    for row in reader:
      if row[pdb_id_index] != current_pdb_id:
        pdb_id = row[pdb_id_index].lower().replace(".pdb", "").replace(".cif", "")
        current_dssr_annotations = pd.read_csv(f"/Users/ibrahims/Documents/Programming/undergrad_reasearch/rna/rna_parsing/{pdb_id}/dssr_annotations.csv")
        current_dssr_annotations["residues_sorted"] = current_dssr_annotations.apply(lambda row: tuple(sorted([str(row["nt1"]).replace("/", ""), str(row["nt2"].replace("/", ""))])), axis = 1)
        current_pdb_id = row[pdb_id_index]
      if row[basepair_index] == "True":
        residue_sorted = tuple(sorted([str(row[nt1_index]), str(row[nt2_index])]))
        try: 
          basepairclassifcation = current_dssr_annotations[current_dssr_annotations["residues_sorted"] == residue_sorted]["LW"].values[0]
        except:
          residue_sorted = tuple(["0" + residue for residue in residue_sorted])
          basepairclassifcation = current_dssr_annotations[current_dssr_annotations["residues_sorted"] == residue_sorted]["LW"].values[0]
      else:
        basepairclassifcation = "None"
      writer.writerow(row + [basepairclassifcation])
  

addclassificationtodata()