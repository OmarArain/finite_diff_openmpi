import numpy as np
import matplotlib.pyplot as plt

outputdir = "output/"
imagedir = "images/"
file = "heat_output_0"
mat_1d = np.fromfile(outputdir+file+".bin", np.float_, -1, "")
dim_size = np.cbrt(mat_1d.size)
mat_3d = np.reshape(mat_1d, (dim_size, dim_size, dim_size))
slice_2d = mat_3d[int(dim_size)/2]
plt.matshow(slice_2d)
plt.savefig(imagedir+file+".png")
