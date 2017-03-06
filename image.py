import numpy as np
import matplotlib.pyplot as plt
import imageio


# 3d files
filenames = []
for i in range(0,2,1):
	outputdir = "output/"
	imagedir = "images/"
	file = "heat_output_3d_"+str(i)

	mat_1d = np.fromfile(outputdir+file+".bin", np.float_, -1, "")
	dim_size = int(np.cbrt(mat_1d.size))
	print (dim_size)
	mat_3d = np.reshape(mat_1d, (dim_size, dim_size, dim_size))
	slice_2d = mat_3d[1]#mat_3d[int(dim_size/2)]
	plt.matshow(slice_2d)
	plt.savefig(imagedir+file+".png")
	filenames.append(imagedir+file+".png")


# images = []
# for filename in filenames:
#     images.append(imageio.imread(filename))
# imageio.mimsave('heat.gif', images)

#2d files
filenames = []
for i in range(0,2, 1):
	outputdir = "output/"
	imagedir = "images/"
	file = "heat_output_2d_"+str(i)

	mat_1d = np.fromfile(outputdir+file+".bin", np.float_, -1, "")
	dim_size = int(np.sqrt(mat_1d.size))
	print (dim_size)
	mat_2d = np.reshape(mat_1d, (dim_size, dim_size))
	#slice_2d = mat_3d[int(dim_size/2)]
	plt.matshow(mat_2d)
	plt.savefig(imagedir+file+".png")
	filenames.append(imagedir+file+".png")


# images = []
# for filename in filenames:
#     images.append(imageio.imread(filename))
# imageio.mimsave('heat.gif', images)
