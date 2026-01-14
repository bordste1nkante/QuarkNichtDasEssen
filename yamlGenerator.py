import yaml

#to store the params of different projects and to make them easily readable by our simulation
projectName= input("Enter yaml file name: ")
h5Name= input("Enter .h5 file name: ")
NumberOfConfigurations = int(input("Enter number of configurations: "))
NumberOfThermalSweeps = int(input("How many calibration steps should be taken: "))
SweepFactor = int(input("How many times should the lattice be updated, before a new configuration is saved: "))
roundingFactor = int(input("After how many global updates should det = 1 be restored: "))
XUpdate = int(input("After how many global updates should the X be reset: "))
xAxis = int(input("Dimension of X-axis: "))
yAxis = int(input("Dimension of Y-axis: "))
zAxis = int(input("Dimension of Z-axis: "))
tAxis = int(input("Dimension of t-axis: "))

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

    "lattice":{
        "x": xAxis,
        "y": yAxis,
        "z":zAxis,
        "t":tAxis
    
    },
    "startConfig": start,
    "updates":{
        "NConfig": NumberOfConfigurations,
        "NSweepsThermal": NumberOfThermalSweeps,
        "Sweep": SweepFactor,
        "Rounding": roundingFactor,
        "XUpdate": XUpdate
    },
    "h5Name":f"{h5Name}.h5"

}


with open(f"{projectName}.yaml", "w") as f:
    yaml.dump(data, f)
