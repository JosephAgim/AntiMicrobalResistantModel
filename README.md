# Project Title

RandomForestClassifier Ability & Limitation Profiling Antibiotic Resistant Protein Structure of Prokaryotes Through Protein Structure Recognition: Antibiotic vs non-Antibiotic Resistant Protein Structures From The Prokaryotic PDB database, A Superviced Machine Learning Approach.

Final project for the Building AI course

## Summary
Antibiotic resistant structural characteristics is believed to have a different atributes than the other proteins existing in the prokaryotic cell line. I tested this theory by building a model that overfitted for the antibiotic resistant protein while being biased for the non-antibiotic resistant proteins(proteins that are not antibiotic resistant).
The code was built through a guide viable at https://biopython.org/. The traning library was downloaded from the https://www.rcsb.org/. Approximatley 80000 PDB files where downloaded in total 70000 was normal proteins in the procaryotic region the other 10000 was antibiotic resistant proteins,  after cleaning by removing the corrupt and duplicate files from the both library there was 8000 PDB files of antibiotic resistant protein structures while having 59000 PDBfies the normal protein structures. Nomalization was used to  by looping the traning of the smaller library 7 times to fit the volume of 59000 since 8000 * 7 = 56000.
The model is overfittet to the traning data. The overall accuracy is 4.53% for protein structures not trained on. Thus the models bias is still towards the normal proteins and not the antibiotic resistant protein structure since the false positives of the normal proteins was 3.40%   

# Contents
Installation
Usage
Data
Contributing
License
Installation

To run the code in this project, follow these steps:

Install python programming language (version 3.14.1) from the python.org website.
Install Visual Studio Code
Install the required Python packages packages by running the following command in the terminal console:

## Requirements
Download the traning library and the model through this link from my google drive https://drive.google.com/drive/folders/1hBu-3AaGMWoH4AqA-biEuj9WN1h6CTtX?usp=sharing
To install the required packages, you can use the provided bash script.

## Installation

1. Clone the repository:

```bash
git clone https://github.com/JosephAgim/AntiMicrobalResistantModel
cd AntiMicrobalResistantModel
bash install_dependencies.sh
```
# Usage
To run the model on your own PDB files you can  use the code below after you change the scripts directory for your own directory containing your files 

```bash
python ModelPlay.py
```

# Contributing
Contributions to this project are welcome. If you would like to contribute, please follow these guidelines:

Fork the repository and create a new branch for your contribution.
Make your changes and submit a pull request.
Ensure your code follows the project's coding style and conventions.
Provide a clear description of your changes and the problem they solve.
For any issues or questions, please use the GitHub issue tracker.

# License
This project is licensed under the MIT License.

# Contact
For any inquiries or feedback, please contact JosephAgim at Josag897@student.liu.se

# Acknowledgements
I would like to acknowledge the invaluable assistance of Uppsala University for providing the computational power through a cluster called uppmax. 

