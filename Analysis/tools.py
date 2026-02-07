import numpy as np

#bootrap method
def bootstrap(iterations, data):
    
    #create from data by drawing with replace x sets (iterations) with len(data) data points
    boot_sample = np.random.choice(data, size=(iterations, len(data)), replace=True)
    #determine mean of each set
    boot_means = np.mean(boot_sample, axis=1)
    #determine overall mean
    boot_mean = np.mean(boot_means)
    #determine standard deviation
    boot_std = np.std(boot_means, ddof=1)

    return boot_mean, boot_std
