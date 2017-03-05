import numpy as np
import matplotlib.pyplot as plt
import imageio

filenames = []
for i in range(0,99):
	outputdir = "output/"
	imagedir = "images/"
	file = "heat_output_"+str(i)

	# mat_1d = np.fromfile(outputdir+file+".bin", np.float_, -1, "")
	# dim_size = np.cbrt(mat_1d.size)
	# mat_3d = np.reshape(mat_1d, (dim_size, dim_size, dim_size))
	# slice_2d = mat_3d[int(dim_size)/2]
	# plt.matshow(slice_2d)
	# plt.savefig(imagedir+file+".png")
	filenames.append(imagedir+file+".png")


images = []
for filename in filenames:
    images.append(imageio.imread(filename))
imageio.mimsave('heat.gif', images)
