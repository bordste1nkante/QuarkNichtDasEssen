import numpy as np
import matplotlib.pyplot as plt
import h5py

#this is local code
import tools


#number of configs
n=2000

#number of bootstrap iterations
B=200


MeanCold = []
StdCold=[]
#read out all data into array
with h5py.File("../h5/ColdPlaq.h5", "r") as f:
    for i in range(n):
        dataset = f[f"Configuration/{i}/plaquette"]
        print(dataset[3])
        mean, std = tools.bootstrap(B, dataset)
        MeanCold.append(mean)
        StdCold.append(std)

MeanHot = []
StdHot=[]
#read out all data into array
with h5py.File("../h5/HotPlaq.h5", "r") as f:
    for i in range(n):
        dataset = f[f"Configuration/{i}"]

        mean, std = tools.bootstrap(B, dataset)
        MeanHot.append(mean)
        StdHot.append(std)


xCold= np.arange(0,len(MeanCold),1)
xHot= np.arange(0,len(MeanHot),1)

plt.errorbar(xCold, MeanCold, yerr = StdCold, label="Cold")
plt.errorbar(xHot, MeanHot, yerr = StdHot, label="Hot")

plt.show()








