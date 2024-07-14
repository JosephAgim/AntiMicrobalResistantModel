import os
import gzip
import numpy as np
from Bio.PDB import PDBParser
import joblib
import warnings
from tqdm import tqdm
from multiprocessing import Pool

warnings.filterwarnings("ignore")

# Convert three-letter residue names to single-letter codes
def three_to_one_letter(residue_name):
    conversion = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    return conversion.get(residue_name, 'X')

# Function to calculate features from a protein structure
def calculate_features(structure):
    hydrophobic_residues = ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y', 'P', 'G']
    hydrophilic_residues = ['R', 'K', 'D', 'E', 'Q', 'N', 'H', 'S', 'T', 'C']
    hydrophobic_count = 0
    hydrophilic_count = 0
    atom_count = 0
    total_residue_count = 0
    distances = []

    for model in structure:
        for chain in model:
            residues = list(chain)
            total_residue_count += len(residues)
            for residue in residues:
                atom_count += len(list(residue.get_atoms()))
                resname = three_to_one_letter(residue.get_resname())
                if resname in hydrophobic_residues:
                    hydrophobic_count += 1
                elif resname in hydrophilic_residues:
                    hydrophilic_count += 1
            for i in range(len(residues) - 2):
                for atom1 in residues[i]:
                    for atom2 in residues[i + 1]:
                        for atom3 in residues[i + 2]:
                            distance1 = atom1 - atom2
                            distance2 = atom2 - atom3
                            distances.append(distance1)
                            distances.append(distance2)

    if distances:
        mean_distance = np.mean(distances)
        std_distance = np.std(distances)
        min_distance = np.min(distances)
        max_distance = np.max(distances)
    else:
        return np.array([0] * 7)  # Fallback for structures with no measurable distances

    hydrophobic_proportion = hydrophobic_count / total_residue_count if total_residue_count else 0
    hydrophilic_proportion = hydrophilic_count / total_residue_count if total_residue_count else 0

    return np.array([mean_distance, std_distance, min_distance, max_distance, hydrophobic_proportion, hydrophilic_proportion, atom_count])

# Main function to process PDB files and classify them using the trained model
def process_pdb_file(pdb_file):
    try:
        with gzip.open(os.path.join(pdb_directory, pdb_file), 'rt') as handle:
            structure = PDBParser().get_structure("new_protein", handle)
            features = calculate_features(structure).reshape(1, -1)
            prediction_proba = clf.predict_proba(features)[0]
            predicted_label = 1 if prediction_proba[1] > 0.5 else 0  # Predicted label based on the higher probability

            # Log output for each file including the predicted label
            print(f"File: {pdb_file}")
            print(f"Predicted probabilities: Non-Resistant: {prediction_proba[0]:.4f}, Resistant: {prediction_proba[1]:.4f}")
            print(f"Predicted label: {'Resistant' if predicted_label == 1 else 'Non-Resistant'}\n")
            
            # Return true if predicted as resistant (label 1), this is useful if you're tracking performance in some way
            return predicted_label == 1
    except Exception as e:
        print(f"Error processing file {pdb_file}: {e}")
        return False

if __name__ == "__main__":
    model_path = "/proj/trained_model.joblib"
    pdb_directory = "/proj/Anti_bio_res_prot"
    clf = joblib.load(model_path)

    pdb_files = [f for f in os.listdir(pdb_directory) if f.endswith('.pdb.gz')]
    with Pool(processes=os.cpu_count()) as pool:
        results = list(tqdm(pool.imap(process_pdb_file, pdb_files), total=len(pdb_files)))

    correct_predictions = sum(results)
    accuracy = correct_predictions / len(results) if results else 0
    print(f"Overall accuracy: {accuracy * 100:.2f}%")
