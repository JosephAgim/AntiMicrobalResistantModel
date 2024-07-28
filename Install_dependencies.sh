
### install_dependencies.sh
#!/bin/bash

# Update package list and install Python3 and pip if not already installed
sudo apt-get update
sudo apt-get install -y python3 python3-pip

# Install Python packages
pip3 install numpy
pip3 install joblib
pip3 install biopython
pip3 install scikit-learn
pip3 install tqdm

echo "All dependencies have been installed."
