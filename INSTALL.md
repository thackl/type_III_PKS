* Create a clean conda environment with Python 3.8

```
conda create -n T3PKS python=3.8
conda activate T3PKS
```

* Intalls dependecies

```
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

  
