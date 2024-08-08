#pymol script for residues visualisation
import pandas as pd
residues = pd.read_csv('../output_data/important_residues.csv', sep = ";")
residues = residues.iloc[:,4:]
from pymol import cmd
import os
#coordinates for selection
x=0
y=0
for item in residues:
    file = "../output_data/AF_pred/residue_importance/" + item + "_importance.pdb"
    cmd.load(file, item)
    #select residues from df
    selector = str()
    for residue in residues[item]:
        if residue != "-":
            selector = selector + " + resi " + str(residue)
    selector = selector[2:]
    selector = "(" + selector + ") and " + item
    selector_name = item + "_residues"
    cmd.select(selector_name, selector)
    
#alignment agaisnt oras
align_list = list(residues.columns)
align_list.remove("ORAS")
for item in align_list:
    cmd.super(item, "ORAS")

#making grid
for item in residues:
    if list(residues.columns).index(item) % 5 ==0:
        x = 0 
        y += -70
        cmd.translate(vector = [x,y,0], selection = item)
    else:
        x += 70
        cmd.translate(vector = [x,y,0], selection = item)
    foo = item + "_foo"
    cmd.pseudoatom(foo, pos = [x+10, y + 20, 0], label = item)
