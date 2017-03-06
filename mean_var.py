import numpy as np

mean_file = "output/heat_output_mean.bin"
var_file = "output/heat_output_std.bin" #wrongly named, they are  actually vars
timesteps = 200
processes = 512
n = 400*400*400/processes
means_1d = np.fromfile(mean_file, np.float_, -1, "")
vars_1d = np.fromfile(var_file, np.float_, -1, "")

means_2d = (np.reshape(means_1d, (processes, timesteps))).T
vars_2d = (np.reshape(vars_1d, (processes, timesteps))).T

pop_means = []
pop_stds = []
for i in range(0,timesteps):
	means_sq = np.power(means_2d[i],2)
	mean = (means_2d[i].sum())/processes
	t = n*(vars_2d[i] + means_sq)
	std = ((t.sum()/(processes*n))-mean**2)**.5
	pop_means.append(mean)
	pop_stds.append(std)


print("MEANS")
pop_means
print("STDS")
pop_stds