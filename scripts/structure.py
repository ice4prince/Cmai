from Bio.PDB import PDBParser
import numpy as np
#read pdb file
def read_pdb(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('antigen',pdb_file)
    return structure

#extract residues and ca atoms
def extract_residues_ca(structure):
    residues = []
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                residues.append(residue)
                ca_atoms.append(residue['CA'])
    return residues, ca_atoms

#compute distance matrix, multiply the matrix by ten
def compute_distance_matrix(ca_atoms):
    distance_matrix = np.zeros((len(ca_atoms),len(ca_atoms)))
    for i in range(len(ca_atoms)):
        for j in range(len(ca_atoms)):
            distance_matrix[i][j] = np.linalg.norm(ca_atoms[i]-ca_atoms[j]) * 10
    return distance_matrix

#compute correlation, convert to vector
def compute_correlation(matrix1, matrix2):
    triu_idx = np.triu_indices(matrix1.shape[0], k=1)
    matrix1_values = matrix1[triu_idx]
    matrix2_values = matrix2[triu_idx]
    correlation = np.corrcoef(matrix1_values, matrix2_values)[0, 1]

    return correlation

#compute rmsd
def compute_rmsd(ca_atoms1, ca_atoms2):
    diff = ca_atoms1 - ca_atoms2
    squared_diff = np.square(diff)
    mean_squared_diff = np.mean(squared_diff)
    rmsd = np.sqrt(mean_squared_diff)
    return rmsd


def analyze_pdb_files(pdb_file1, pdb_file2):
    antigen1 = read_pdb(pdb_file1)
    antigen2 = read_pdb(pdb_file2)

    a1_residues, a1_ca_atoms = extract_residues_ca(antigen1)
    a2_residues, a2_ca_atoms = extract_residues_ca(antigen2)

    if len(a1_residues) != len(a2_residues):
        raise ValueError('The two antigens have different number of residues.')

    a1_distance_matrix = compute_distance_matrix(a1_ca_atoms)
    a2_distance_matrix = compute_distance_matrix(a2_ca_atoms)

#    rmsd = compute_rmsd(a1_distance_matrix, a2_distance_matrix)
    correlation = compute_correlation(a1_distance_matrix, a2_distance_matrix)

#    return rmsd, correlation
    return correlation

print(analyze_pdb_files('ranked_1.pdb','ranked_1.pdb'))



