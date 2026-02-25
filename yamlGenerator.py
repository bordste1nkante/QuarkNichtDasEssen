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
xAxis = int(input("Dimension of X-axis: "))
yAxis = int(input("Dimension of Y-axis: "))
zAxis = int(input("Dimension of Z-axis: "))
tAxis = int(input("Dimension of t-axis: "))


NumberOfConfigurations = int(input("Enter number of configurations: "))
NumberOfThermalSweeps = int(input("How many calibration steps should be taken: "))
SweepFactor = int(input("How many times should the lattice be updated, before a new configuration is saved: "))
XUpdate = int(input("After how many global updates should the X be reset: "))
Multi = int(input("How many multihits: "))
over = int(input("After how many steps a overrelaxation step: "))
epsilon = float(input("What epsilon value: "))

beta = float(input("Beta: "))
testrun = bool(input("Is this a test or debugging run, 1 yes 0 no: "))

#seed1 = random.uniform(0,2)
seed2 = random.uniform(100,200)
#seed3 = random.uniform(0,4)
seed4 = random.uniform(250,350)
#seed5 = random.uniform(0,6)

#seed3_1 = random.uniform(7,10 )
#seed3_2=random.uniform(7, 11)
#seed3_3 = random.uniform(7, 12)
#seed3_4 =random.uniform(7, 13)
#seed3_5 = random.uniform(7, 14)


seed3_6  = random.uniform(400, 500)
#seed3_7  = random.uniform(10, 22)
#seed3_8  = random.uniform(10, 23)
seed3_9  = random.uniform(70, 300)
#seed3_10= random.uniform(10, 26)
#seed3_11= random.uniform(10, 27)
seed3_12= random.uniform(173, 350)
#seed3_13= random.uniform(10, 29)
#seed3_14= random.uniform(10, 24)

#seedOver = random.uniform (100,102)


if (testrun==True):
    #seed1 = 69
    seed2 = 42
    #seed2_1 = 420
    #seed2_2 = 4200

    #seed3 = 67
    #seed3_1 = 1
    #seed3_2=2
    #seed3_3 = 3
    #seed3_4 = 4
    #seed3_5 = 5

    seed3_6  = 10
    #seed3_7  = 11
    #seed3_8  = 12
    seed3_9  = 20
    #seed3_10= 21
    #seed3_11= 22
    seed3_12= 30
    #seed3_13= 31
    #seed3_14= 32

    seed4 = 666
    #seed5 = 1949

    #seedOver = 2026
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
    
    },
    "startConfig": start,
    "epsilon": epsilon,
    "updates":{
        "NConfig": NumberOfConfigurations,
        "NSweepsThermal": NumberOfThermalSweeps,
        "Sweep": SweepFactor,
        "XUpdate": XUpdate,
        "Multi": Multi,
        "overrelaxation": over
    },
    "seeds":{
        #"distEpsilon": seed1,
        "dist": {"rng1": seed2,
                # "rng2": seed2_1,
                # "rng3": seed2_2,
                # "hotNumb1": seed3,                
                # "hotNumb2": seed3_1,
                # "hotNumb3":seed3_2,
                # "hotNumb1Extra" : seed3_3,               
                # "hotNumb2Extra" :seed3_4,               
                # "hotNumb3Extra"  :seed3_5,

                 "hotNumbT1": seed3_6,
                 #"hotNumbT2": seed3_7,
                 #"hotNumbT3": seed3_8,
                 "hotNumbS1": seed3_9,
                 #"hotNumbS2": seed3_10,
                 #"hotNumbS3": seed3_11,
                 "hotNumbR1": seed3_12,
                 #"hotNumbR2": seed3_13,
                 #"hotNumbR3": seed3_14
                 },
        "indexDist": seed4,
        #"uniformAcceptReject": seed5,
        #"overrelaxation": seedOver
    },
    "h5Name":f"{projectName}.h5"

}


with open(f"{projectName}.yaml", "w") as f:
    yaml.dump(data, f)
