import yaml
import numpy as np
from numpy import random


#to store the params of different projects and to make them easily readable by our simulation
projectName= input("Enter yaml file name: ")
#h5Name= input("Enter .h5 file name: ")
observable = int(input("Enter observable code: "))
startPoint = np.zeros(4, dtype=np.int64)
endPoint = np.zeros(4, dtype=np.int64)
x= ["x", "y", "z", "t"]
if(observable == 1 or observable == 2):
    for i in range(0,4):
        startPoint[i]= int(input(f"Enter start position of {x[i]}-axis: "))
        endPoint[i] =  int(input(f"Enter end position of {x[i]}-axis: "))



NumberOfConfigurations = int(input("Enter number of configurations: "))
NumberOfThermalSweeps = int(input("How many calibration steps should be taken: "))
SweepFactor = int(input("How many times should the lattice be updated, before a new configuration is saved: "))
roundingFactor = int(input("After how many global updates should det = 1 be restored: "))
XUpdate = int(input("After how many global updates should the X be reset: "))
xAxis = int(input("Dimension of X-axis: "))
yAxis = int(input("Dimension of Y-axis: "))
zAxis = int(input("Dimension of Z-axis: "))
tAxis = int(input("Dimension of t-axis: "))
beta = float(input("Coupling value: "))
a = float(input("Lattice spacing: "))
testrun = bool(input("Is this a test or debugging run, 1 yes 0 no: "))

seed1 = random.uniform(0,2)
seed2 = random.uniform(0,3)
seed3 = random.uniform(0,4)
seed4 = random.uniform(0,5)
seed5 = random.uniform(0,6)

if (testrun==True):
    seed1 = 69
    seed2 = 42
    seed3 = 67
    seed4 = 666
    seed5 = 1949

loop = True
start = True
while(loop):
    startTemp = int(input("Should there be a cold or a hot start, enter 1 or 0: "))
    if(startTemp==0):
        start = False
        loop = False
    if(startTemp ==1):
        loop = False



data= {

    "observable": observable,
    "positions":{
        "startPoint": startPoint.tolist(),
        "endPoint": endPoint.tolist()
    },
    "lattice":{
        "x": xAxis,
        "y": yAxis,
        "z":zAxis,
        "t":tAxis,
        "beta": beta,
        "lattice spacing": a
    
    },
    "startConfig": start,
    "updates":{
        "NConfig": NumberOfConfigurations,
        "NSweepsThermal": NumberOfThermalSweeps,
        "Sweep": SweepFactor,
        "Rounding": roundingFactor,
        "XUpdate": XUpdate
    },
    "seeds":{
        "distEpsilon": seed1,
        "dist": seed2,
        "hotDist": seed3,
        "indexDist": seed4,
        "uniformAcceptReject": seed5
    },
    "h5Name":f"{projectName}.h5"

}


with open(f"{projectName}.yaml", "w") as f:
    yaml.dump(data, f)
