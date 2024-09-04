# type_III_PKS

### predictive model
#create a new conda environment 
conda create -n "ml" python=3.8
conda activate ml

#install jupyter notebook
conda install anaconda::jupyter

#cd to the folder with the script and run it
jupyter notebook Predictive_ML.ipynb

#install the necessary dependencies
conda install anaconda::h5py
pip install rdkit
conda install conda-forge::matplotlib
conda install anaconda::scikit-learn

#added "import random" in cell 6
#rdkit needs to be installed with pip, conda doesn't work
#cell 8: "training_all_score = score_all_points(score_all)"

### descriptive model
#install biopython in jupyter notebook first
#add "T3PKS_alignment_clustal" to the input folder
