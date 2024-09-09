# Description

The repository contains datasets and scripts to reproduce results and figures from the paper _Nika Sokolova, Stepan S. Denisov, Thomas Hackl, and Kristina Haslinger_
_**"Unravelling the functional diversity of type III polyketide synthases in fungi "**_


The code was implemented on MacOS with the following crucial packages:

*Python 3.8.
*h5py 3.7.0
*biopython 1.79
*jupyterlab 4.0.10
*numpy 1.24.4
*pandas 1.5.3
*peptides 0.3.2
*rdkit 2023.3.3
*scikit-learn 0.24.2

Packages could be installed using `pip` and `conda`
```
conda create -n T3PKS python=3.8
conda activate T3PKS
pip install jupyterlab
pip install rdkit
pip install h5py
conda install anaconda::scikit-learn
conda install conda-forge::matplotlib
pip install peptides
pip install biopython
```

* Download and unzip type_III_PKS folder and run Jupyter

```
git clone https://github.com/denisovss/type_III_PKS.git
cd type_III_PKS
jupyter notebook
```

