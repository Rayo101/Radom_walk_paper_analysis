'''
The following code converts the *.root ouput of the GATE simuation and extracts the pair of locations of gamma 
ray interactions in the crystals as well of the time of both of these interactions and writes them to a *.npy 
file that contains the above information.  This code was based off the listmode_ROOT2CSV.py script developed
for the G4PEPT Geant4 simulation.  

Author:            Rayhaan Perin 
Last Update:       14-08-2021
Based On:          listmode_ROOT2CSV.py from G4PEPT
Required Packages: uproot, numpy

To install packages packages use pip3 install package or pip install package,
alternatively if one uses conda you can run conda install package
'''

####################
######PACKAGES######
####################

import uproot
import numpy as np
import os 
from natsort import os_sorted
from tqdm import tqdm

##################
######INPUTS######  
##################

INPUT_PATH = "/home/rayhaan/TimingInvestigations/TimingInvestigations/GateData/root/"
COINCIDENCES_TREE= "Coincidences"
SINGLES_TREE = "Singles"
OUT_PATH = "/home/rayhaan/TimingInvestigations/TimingInvestigations/GateData/raw/"

#####################
######FUNCTIONS######
#####################

def moveToFace(x,y, distance):
    x = x + distance*np.cos(np.arctan2(y,x))
    y = y + distance*np.sin(np.arctan2(y,x))
    return x, y 

def folderIfNotExist(path):
    if not os.path.isdir(path):
        os.makedirs(path)

#####################
######ROOT2CSV#######
#####################

folderFiles = np.array(os_sorted(os.listdir(INPUT_PATH))).astype(np.str0)
print("Looping through Files: \n")
for i in tqdm(range(len(folderFiles))):

    MyFile = uproot.open(INPUT_PATH + folderFiles[i])
    CoincidenceTree = MyFile[COINCIDENCES_TREE]

    headers = ["xA","yA","zA","xB","yB","zB","tA","tB"]#,"xg","yg","zg","displacement","PositronInitialEnergy"]

    # Assign a variable to each leaf in the reconData tree
    xA_leaf = "globalPosX1"
    yA_leaf = "globalPosY1"
    zA_leaf = "globalPosZ1"
    xB_leaf = "globalPosX2"
    yB_leaf = "globalPosY2"
    zB_leaf = "globalPosZ2"
    tB_leaf = "time2"
    tA_leaf = "time1"

    # Form an array for each leaf
    x1 = CoincidenceTree[xA_leaf].array()
    y1 = CoincidenceTree[yA_leaf].array()
    z1 = CoincidenceTree[zA_leaf].array()
    x2 = CoincidenceTree[xB_leaf].array()
    y2 = CoincidenceTree[yB_leaf].array()
    z2 = CoincidenceTree[zB_leaf].array()
    t2 = CoincidenceTree[tB_leaf].array()
    t1 = CoincidenceTree[tA_leaf].array()

    # Move the simulated interaction positions to the face of the crystal.  The interation position is moved 5 mm radially inwards.  
    x1, y1 = moveToFace(x1,y1, -5)
    x2, y2 = moveToFace(x2,y2, -5)

    # Write the formatted data  
    saveArr = np.array([x1, y1, z1, x2, y2, z2, t1, t2], dtype = np.float64).T
    np.save(arr = saveArr, file = "{}".format(OUT_PATH + os.path.splitext(folderFiles[i])[0]))

