"""
======================================================================
Molecular Module Demo
======================================================================

A small example to show how to use the various molecular
representations present in the `molecular` module to visualize
proteins. This example also shows how to parse PDB files to obtain
atomic info essential for constructing the representations.

Importing necessary modules
"""

import urllib
import os
from fury import window, actor, ui, molecular as mol
import numpy as np

###############################################################################
# Downloading the PDB file of the protein to be rendered.
# User can change the pdb_code depending on which protein they want to
# visualize.
pdb_code = '1crn'
downloadurl = "https://files.rcsb.org/download/"
pdbfn = pdb_code + ".pdb"
flag = 0
if not os.path.isfile(pdbfn):
    flag = 1
    url = downloadurl + pdbfn
    outfnm = os.path.join(pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
    except Exception:
        print("Error in downloading the file!")

###############################################################################
# creating a PeriodicTable() object to obtain atomic numbers from names of
# elements
table = mol.PeriodicTable()

###############################################################################
# Creating empty lists which will be filled with atomic information as we
# parse the pdb file.
NumberOfAtoms = 0

points = []
elements = []
atom_names = []
model = []
sheets = []
helix = []
residue_seq = []
chain = []
is_hetatm = []
current_model_number = 1

###############################################################################
# Parsing the pdb file for information about coordinates and atoms

pdbfile = open(pdbfn, 'r')
pdb_lines = pdbfile.readlines()
for line in pdb_lines:
    line = line.split()
    try:
        if line[0] == 'ATOM' or line[0] == 'HETATM':
            if line[-1] != 'H':
                coorX, coorY, coorZ = float(line[6]), float(line[7]), \
                                      float(line[8])
                resi = line[5]
                current_chain = ord(line[4])
                points += [[coorX, coorY, coorZ]]
                residue_seq += [resi]
                chain += [current_chain]
                elements += [table.atomic_number(line[-1])]
                atom_names += [line[2]]
                model += [current_model_number]
                NumberOfAtoms += 1
                if(line[0] == 'HETATM'):
                    is_hetatm += [True]
                else:
                    is_hetatm += [False]
        if line[0] == 'SHEET':
            start_chain = ord(line[5])
            start_resi = int(line[6])
            end_chain = ord(line[8])
            end_resi = int(line[9])
            r = [start_chain, start_resi, end_chain, end_resi]
            sheets += [r]
        if line[0] == 'HELIX':
            start_chain = ord(line[4])
            start_resi = int(line[5])
            end_chain = ord(line[7])
            end_resi = int(line[8])
            r = [start_chain, start_resi, end_chain, end_resi]
            helix += [r]
    except Exception:
        continue

points = np.array(points)
residue_seq = np.array(residue_seq, dtype=int)
chain = np.array(chain)
elements = np.array(elements)
atom_names = np.array(atom_names)
model = np.array(model)
sheets = np.array(sheets)
helix = np.array(helix)
is_hetatm = np.array(is_hetatm)


###############################################################################
# Helper function to make the visuals look good by manipulating lighting.
def make_aesthetic(molecule_rep):
    molecule_rep.GetProperty().SetAmbient(0.2)
    molecule_rep.GetProperty().SetDiffuse(1)
    molecule_rep.GetProperty().SetSpecular(1)
    molecule_rep.GetProperty().SetSpecularPower(100.0)


###############################################################################
# Doing 3 things here -
# 1. Creating the molecule object.
# 2. Computing the bonding information for the molecule.
# 3. Generating and adding molecular representations to the scene.

molecule = mol.Molecule(elements, points, atom_names, model,
                        residue_seq, chain, sheets, helix, is_hetatm)
# print(np.unique(elements))
mol.compute_bonding(molecule)

# bounding box
b_box = mol.bounding_box(molecule, colors=(0, 0.8, 1), linewidth=0.4)

# stick representation
stick_rep = mol.stick(molecule, bond_thickness=0.2)
make_aesthetic(stick_rep)

# ribbon representation
ribbon_rep = mol.ribbon(molecule)
make_aesthetic(ribbon_rep)

# ball and stick representation
ball_stick_rep = mol.ball_stick(molecule, atom_scale_factor=0.3,
                                bond_thickness=0.2)
make_aesthetic(ball_stick_rep)

# sphere representation
vdw_sphere_rep = mol.sphere_cpk(molecule)
make_aesthetic(vdw_sphere_rep)

###############################################################################
# Dimensions of the output screen
screen_x_dim = 800
screen_y_dim = 300
dims = (screen_x_dim, screen_y_dim)

###############################################################################
# creating a ShowManager object
showm = window.ShowManager(size=dims, title=pdb_code)


tb = ui.TextBlock2D(text=pdb_code.upper(), position=(screen_x_dim/2-40,
                    screen_y_dim/30), font_size=20, color=(1, 1, 1))
tb.actor.GetTextProperty().SetFontFamilyToCourier()

###############################################################################
# Creating the textblocks and molecular representations.

reps = [stick_rep, ribbon_rep, ball_stick_rep, vdw_sphere_rep]
rep = ['Stick', 'Ribbon', 'Ball & Stick', 'Sphere']
text = []
for i in range(4):
    t_actor = actor.label(rep[i], pos=(0, 0, 0),
                          scale=(2, 2, 2))
    t_actor.SetCamera(showm.scene.camera())
    text.append(t_actor)

###############################################################################
# Grid for the representations
grid_ui = ui.GridUI(actors=reps, captions=text,
                    caption_offset=(-6, 0, 20), dim=(1, 4),
                    cell_padding=2,
                    aspect_ratio=1,
                    rotation_axis=(0, 1, 0))
showm.scene.add(grid_ui)

showm.scene.reset_clipping_range()
showm.scene.add(tb)

showm.scene.set_camera((39.169, 261.401, 84.727),
                       (51.889, 20.541, 12.977),
                       (0.0, 0.28, -0.96))
showm.scene.zoom(2.5)

###############################################################################
# Delete the PDB file.
flag = 0
if flag:
    os.remove(pdbfn)


interactive = True
if interactive:
    showm.start()

###############################################################################
# to save a snapshot of the image
window.record(showm.scene, size=dims, out_path='images/grid_'+pdb_code+'.png')
