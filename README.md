# type_III_PKS

# create a new conda environment 
conda create -n "ml" python=3.8
conda activate ml

# install jupyter notebook
conda install anaconda::jupyter

# cd to the folder with the script and run it
jupyter notebook Predictive_ML.ipynb

# install the necessary dependencies
conda install anaconda::h5py
conda install conda-forge::rdkit
conda install conda-forge::matplotlib
conda install anaconda::scikit-learn
