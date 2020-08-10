import numpy
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot


# Create a new plot
figure = pyplot.figure()
ax = mplot3d.Axes3D(figure)
ax.set_xlim(-50,50)
ax.set_ylim(-50,50)
ax.set_zlim(-50,50)

#Function to calculate centre of each triangle.

# Load the STL files and add the vectors to the plot
your_mesh = mesh.Mesh.from_file(r'C:\Users\User\Desktop\Slicer 4.10.2\NEW STUFF\UROP extension\Skull_rotated_hardened_top_decimated.stl')
# ax.add_collection3d(mplot3d.art3d.Poly3DCollection(your_mesh.vectors))

#Invert all normal vectors
ngtv_mesh = your_mesh.normals * -1

#Calculates the midpoint of each triangle and plots the normal starting from the midpoint
for i in range(0,len(your_mesh.normals)):
    ave_coord = []
    ave_coord.append((your_mesh.vectors[i,0][0] + your_mesh.vectors[i,1][0] + your_mesh.vectors[i,2][0]) / 3)
    ave_coord.append((your_mesh.vectors[i,0][1] + your_mesh.vectors[i,1][1] + your_mesh.vectors[i,2][1]) / 3)
    ave_coord.append((your_mesh.vectors[i,0][2] + your_mesh.vectors[i,1][2] + your_mesh.vectors[i,2][2]) / 3)
    ax.quiver(ave_coord[0], ave_coord[1], ave_coord[2], ngtv_mesh[i,0], ngtv_mesh[i,1], ngtv_mesh[i,2])

# Auto scale to the mesh size
#scale = your_mesh.points.flatten('C')


pyplot.show()

#print(your_mesh.vectors.shape)
