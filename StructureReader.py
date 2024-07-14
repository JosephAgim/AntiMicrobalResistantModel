import os
import gzip
import numpy as np
import warnings
import joblib
import time
from Bio.PDB import PDBParser
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from Bio import BiopythonWarning
from multiprocessing import Pool
from tqdm import tqdm

warnings.simplefilter('ignore', BiopythonWarning)

# Add a start time
start_time = time.time()

# Convert three-letter amino acid residue names to single-letter amino acid residue codes since PDB files has the three-letter names
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

 
    for model in structure: # The for-loop is iterating the model through the protein structure
        for chain in model: # Iterating the chain in model
            residues = list(chain) # Making residues variable follow the protein chain list 
            total_residue_count += len(residues) # Thus the lenght of the residues would be the total aminoacid residue length of the protein
            for residue in residues: # Iterating residue through residues that is the chian list
                atom_count += len(list(residue.get_atoms())) # The length of atoms per residue in each list total amount of atoms
                resname = three_to_one_letter(residue.get_resname()) # calls function to  convert the 3 letter residue name to one letter code
                if resname in hydrophobic_residues: # Finds the amount of hydrophobic one letter code in each residue
                    hydrophobic_count += 1 # counts the hydrophobic residues
                elif resname in hydrophilic_residues: # Finds the amount of hydrophilic one letter code in each residue
                    hydrophilic_count += 1 # counts the hydrophilic residues
            for i in range(len(residues) - 2): # Iterating i in residue chain of the protein to have process go through all the amino acids of the protein, added -2 since we start on 2
                for atom1 in residues[i]:  # Finding the first atom of each residue(amino acid) 
                    for atom2 in residues[i+1]: # Finding the secomnd atom of each residue(amino acid)
                        for atom3 in residues[i+2]: # Finding the third atom of each residue(amino acid)
                            distance1 = atom1 - atom2 # Calculating the distance between atom1 and atom2
                            distance2 = atom2 - atom3 # Calculating the distance between atom2 and atom3
                            distances.append(distance1) # Filling the distances list with the first distance 
                            distances.append(distance2) # Filling the distance list with  the second distance 


    # Calculate summary statistics for distances
    if distances: #Since the distance list is filled we can let the numpy np module work the mean, std, min, and max value of the list 
        mean_distance = np.mean(distances)
        std_distance = np.std(distances)
        min_distance = np.min(distances)
        max_distance = np.max(distances)
    else: # if ther is no distance in the list
        print(f"No distances found for protein {structure.id}") # Print
        mean_distance = std_distance = min_distance = max_distance = 0 #and make al values 0

    hydrophobic_proportion = hydrophobic_count / total_residue_count if total_residue_count else 0 # the hydrophobic proportion of the residues are calulated 
    hydrophilic_proportion = hydrophilic_count / total_residue_count if total_residue_count else 0 # the hydrophilic proportion of the residues are calulated

    # Return feature array including atom count
    return np.array([mean_distance, std_distance, min_distance, max_distance, hydrophobic_proportion, hydrophilic_proportion, atom_count]) # Save all values in numpy array

# Function to process a single file
def process_file(args, parser): 
    dirName, fname, label = args # variables
    print(f'\tFound PDB file: {fname}') # print the protein name when file is located in the directory
    try:
        with gzip.open(os.path.join(dirName, fname), 'rt') as handle: # Opens compressed PDB File and reads in text 'rt' 
            if handle.read(1):  # Check if file is not empty
                handle.seek(0)  # Reset file pointer to beginning
                structure = parser.get_structure(f"{fname}", handle) # gets structure of the PDB file 
                # Calculate features for this protein
                protein_features = calculate_features(structure) # Using calculate_feature() function on the structure from the PDB file 
                return protein_features, label # returning protein freatures and a label para for args
    except Exception as e:  # Catch all exceptions
        print(f"Error processing file {fname}: {e}")
        return None

def process_file_with_parser(args_parser): # for processing file with parser  
    args, parser = args_parser # variables args for label and parser for args_parser 
    return process_file(args, parser) # return process_file() function while placing parameter args and parser where parser is args_parser 

# Initialize array to store features and labels
features = [] 
labels = []

# Function to process a directory                   
def process_directory(rootDir, label, parser):
    for dirName, _, fileList in os.walk(rootDir): # Iterating three variabels through the os.wals in directory 
        print(f'Found directory: {dirName}') # print directory name when found
        args = [(dirName, fname, label) for fname in fileList if fname.endswith('.pdb.gz')] # adding three variables = args while itterating fname in fileList if the PDB file ends with .pdb.gz for compressed file 
        args_parser = [(arg, parser) for arg in args]  # Pair arg with parser while iterating arg in args pairing args with parser 
        with Pool() as pool: # using mutiprocessing
            for result in tqdm(pool.imap(process_file_with_parser, args_parser), total=len(args), desc="Processing files"): 
                if result and result[0] is not None:  # Check if features are not None
                    features.append(result[0])
                    labels.append(result[1])




# Training Library Directory
Anti_rootDir = "/proj/Anti_bio_res_prot"
non_Anti_rootDir = "/proj/non_Anti_bio_res_prot"


def main():
    global features, labels # makign the empry variaable array global
    try:
        # Add a start time
        start_time = time.time()


        # Create a PDBParser object
        parser = PDBParser()

        # Define the number of repetitions for the antibiotic-resistant directory in order to  normalize the valume of PDB files per label
        repetitions = 6 

        # Process the antibiotic-resistant directory multiple times 
        for _ in range(repetitions):
            process_directory(Anti_rootDir, 1, parser)  # Antibiotic resistant dir, label, amd parser


        process_directory(non_Anti_rootDir, 0, parser)  # Non antibiotic resistant dir, label, and parser

        # Convert features and labels to NumPy arrays
        features = np.array(features)
        labels = np.array(labels)

        # Replace NaN values with the mean of the column
        features = np.where(np.isnan(features), np.ma.array(features, mask=np.isnan(features)).mean(axis=0), features)



        # Check that the features and labels lists have the same length
        if len(features) != len(labels):
            print("Error: features and labels lists have different lengths")

        # Check for missing values in the features
        if np.any(np.isnan(features)):
            print("Error: features array contains missing values")

        # Check that the labels are within the range [0, 1]
        if not all(0 <= label <= 1 for label in labels):
            print("Error: labels are not within the range [0, 1]")

        # Check that all features have the expected shape
        for i, feature in enumerate(features):
            if feature.shape != (7,):  # Expected Shape is 7
                print(f"Error: feature at index {i} has incorrect shape {feature.shape}")


        # Split the data into a 80% training set and a 20% test set. Sets random_state = 42 for reproducable result 
        features_train, features_test, labels_train, labels_test = train_test_split(features, labels, test_size=0.2, random_state=42)

        # Train a random forest classifier on the training data iwth the highest number of estimators 100 
        clf = RandomForestClassifier(n_estimators=100, random_state=42)
        clf.fit(features_train, labels_train)

        # Use the trained classifier to make predictions on the test data
        labels_pred = clf.predict(features_test)
        labels_proba = clf.predict_proba(features_test)[0]
        # Calculate the accuracy of the predictions
        accuracy = accuracy_score(labels_test, labels_pred)
        print(f"Accuracy: {accuracy}")


        # Use the trained classifier to make predictions on the test data
        labels_pred = clf.predict(features_test)

        # Calculate the accuracy of the predictions
        accuracy = accuracy_score(labels_test, labels_pred)
        print(f"Accuracy: {accuracy}")

        # Import necessary libraries and functions
        from sklearn.metrics import (confusion_matrix, precision_score, recall_score, f1_score, roc_auc_score,
                             matthews_corrcoef, mean_squared_error, r2_score,
                             balanced_accuracy_score, log_loss, hamming_loss, jaccard_score,
                             mean_absolute_error, median_absolute_error, explained_variance_score)

        # Print the confusion matrix
        print("Confusion Matrix:")
        print(confusion_matrix(labels_test, labels_pred))

        # Print the precision, recall, and F1 score for classification
        print("Precision:", precision_score(labels_test, labels_pred))
        print("Recall:", recall_score(labels_test, labels_pred))
        print("F1 Score:", f1_score(labels_test, labels_pred))

        # Print the ROC AUC
        # Ensure labels_proba is the probability outputs of your classifier for multiclass handling
        print("AUC:", roc_auc_score(labels_test, labels_pred))

        # Print the Matthews correlation coefficient
        print("MCC:", matthews_corrcoef(labels_test, labels_pred))

        # Print classification metrics
        print("Balanced Accuracy:", balanced_accuracy_score(labels_test, labels_pred))
        print("Jaccard Score:", jaccard_score(labels_test, labels_pred))
        print("Hamming Loss:", hamming_loss(labels_test, labels_pred))
        print("Log Loss:", log_loss(labels_test, labels_pred))

        # Print regression metrics
        print("Mean Squared Error:", mean_squared_error(labels_test, labels_pred))
        print("Mean Absolute Error:", mean_absolute_error(labels_test, labels_pred))
        print("Median Absolute Error:", median_absolute_error(labels_test, labels_pred))
        print("R2 Score:", r2_score(labels_test, labels_pred))
        print("Explained Variance Score:", explained_variance_score(labels_test, labels_pred))
        # Save the trained model to a file
        joblib.dump(clf, "/proj/trained_model.joblib")

        # record End time
        end_time = time.time()

        # Calculate and print time
        time_taken = end_time - start_time
        print(f"Time taken : {time_taken} ")



    finally:
        print("Done")

if __name__ == '__main__':
    main()