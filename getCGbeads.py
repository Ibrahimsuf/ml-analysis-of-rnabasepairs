import pandas as pd
import Bio
from Bio.PDB import PDBParser, MMCIFParser
import os
from preprocess_data import Preprocessor
from Bio.PDB import Atom, Residue


G_atoms = ['P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'O3\'', 'C2\'', 'O2\'', 'C1\'', 'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4']
A_atoms = ['P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', "O4'", 'C3\'', 'O3\'', 'C2\'', 'O2\'', 'C1\'', 'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4']
C_atoms = ['P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'O3\'', 'C2\'', 'O2\'', 'C1\'', 'N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6']
U_atoms = ['P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'O3\'', 'C2\'', 'O2\'', 'C1\'', 'N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6']

atoms_needed = {"G": G_atoms, "A": A_atoms, "C": C_atoms, "U": U_atoms}
  
def getAllBeadLocations(residue_pair_distances_path, pdbs_dir, output_path):
  residue_pair_distances = pd.read_csv(residue_pair_distances_path)
  beads = []
  missing_residues = []
  for pdb_id in residue_pair_distances["pdb_id"].unique():
    print(pdb_id)
    parser = PDBParser() if pdb_id.endswith(".pdb") else MMCIFParser()
    path = os.path.join(pdbs_dir, pdb_id)
    model = parser.get_structure(pdb_id, path)[0]
    for chain in model:
      for residue in chain:
        resname = Preprocessor.replaceHetatms(residue)
        if not resname in ["A", "U", "G", "C"]:
          continue
        atoms_in_residue = set([atom.get_id() for atom in residue.get_atoms()])

        atoms_needed_for_residue = atoms_needed[resname]
        atoms_needed_for_residue = set(atoms_needed_for_residue)

        if not atoms_needed_for_residue.issubset(atoms_in_residue):
          missing_residues.append((pdb_id, chain.id, residue.id, resname))
          continue

        CG_beads = get_CG_Beads(residue, resname)
        CG_beads["pdb_id"] = pdb_id
        CG_beads["chain"] = chain.id
        CG_beads["residue"] = residue.id
        CG_beads["resname"] = resname
        beads.append(CG_beads)

  pd.DataFrame(missing_residues, columns=["pdb_id", "chain", "residue", "resname"]).to_csv("missing_residues.csv")
  pd.DataFrame(beads).to_csv(output_path)

def get_geometric_center(residue, atoms):
  CGbead = Residue.Residue(" ", "bead", " ")

  for atom in atoms:
    CGbead.add(residue[atom])

  return tuple(CGbead.center_of_mass(geometric=True))

def get_CG_Beads(residue, resname):
  P_coords = get_geometric_center(residue, ["P", "OP1", "OP2", "O5'", "O3'"])
  S_coords = get_geometric_center(residue, ["C5'", "C4'", "O4'", "C3'", "C2'", "O2'", "C1'"])

  if resname in ["A", "G"]:
    R1_coords = get_geometric_center(residue, ["N9", "C8", "N7", "C5", "C4", "N3"])

  if resname == "A":
    A1_coords = get_geometric_center(residue, ["C6", 'N6'])
    A2_coords = get_geometric_center(residue, ["C2", "N1"])

  if resname == "G":
    G1_coords = get_geometric_center(residue, ["C6", "O6", "N1"])
    G2_coords = get_geometric_center(residue, ["C2", "N2"])

  if resname in ["C", "U"]:
    Y1_coords = get_geometric_center(residue, ["N1", "C5", "C6"])
    Y2_coords = get_geometric_center(residue, ["C2", "O2"])
  if resname == "C":
    C1_coords = get_geometric_center(residue, ["N3", "C4", "N4"])
  if resname == "U":
    U1_coords = get_geometric_center(residue, ["N3", "C4", "O4"])

  CG_beads = {}
  variables = {"P_coords": "P_coords", "S_coords": "S_coords", "R1_coords": "R1_coords", "A1_coords": "A1_coords", "A2_coords": "A2_coords", "G1_coords": "G1_coords", "G2_coords": "G2_coords", "Y1_coords": "Y1_coords", "Y2_coords": "Y2_coords", "C1_coords": "C1_coords", "U1_coords": "U1_coords"}

  for key, value in variables.items():
      if value in locals():
          CG_beads[key] = locals().get(value, globals().get(value))

  return CG_beads

def main():
  getAllBeadLocations("3_5/residue_pair_distances.csv", "/Users/ibrahims/Documents/Programming/undergrad_reasearch/rna/rna_annotations_parser/pdbs/3_5/cif_files", "3_5/beads.csv")


if __name__ == "__main__":
  main()