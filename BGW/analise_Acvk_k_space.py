
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

# exciton to be analized
iexc = int(sys.argv[1])
exciton_file = sys.argv[2]
plot_periodic_images = False
put_everything_1BZ = False

if len(sys.argv) == 4:
    if sys.argv[3] == 'True':
        print('Plotting periodic images of Acvk beyond the 1st BZ!')
        plot_periodic_images = True
    else:
        print('Just plotting Acvk on the 1st BZ!')

def get_exciton_info(exciton_file, iexc):

    """    
    Return the exciton energy and the eigenvec coefficients Acvk

    Assuming calculations with TD approximation
    Info about file at: http://manual.berkeleygw.org/3.0/eigenvectors_h5_spec/
    Also, just working for excitons with Q = 0 and one spin

    TODO -> for now calculting exciton info for exciton with index iexc
    but later, make it calculate for and set of exciton indexes
    
    Parameters:
    exciton_file = exciton file name (string). ex: eigenvecs.h5
    iexc = Exciton index to be read
    
    Returns:
    Acvk = Exciton wavefunc coefficients. array Akcv[ik, ic, iv] with complex values
    Omega = Exciton energy (BSE eigenvalue) in eV (float)
    """

    print('Reading exciton info from file', exciton_file)

    f_hdf5 = h5py.File(exciton_file, 'r')
    
    flavor_calc = f_hdf5['/exciton_header/flavor'][()]
    eigenvecs   = f_hdf5['exciton_data/eigenvectors'][()]          # (nQ, Nevecs, nk, nc, nv, ns, real or imag part)
    eigenvals   = f_hdf5['exciton_data/eigenvalues'][()] 
    
    if flavor_calc == 2:
        print('Flavor in BGW: complex')
        Akcv = eigenvecs[0,iexc-1,:,:,:,0,0] + 1.0j*eigenvecs[0,iexc-1,:,:,:,0,1]
    else:
        print('Flavor in BGW: real')
        Akcv = eigenvecs[0,iexc-1,:,:,:,0,0]

    # K points in the fine grid
    Kpoints_bse = f_hdf5['/exciton_header/kpoints/kpts'][()] 
    Nkpoints = f_hdf5['/exciton_header/kpoints/nk'][()]

    return Akcv, Kpoints_bse, Nkpoints


# getting exciton info
Akcv, Kpoints_bse, Nkpoints = get_exciton_info(exciton_file, iexc)

# assuming a uniform grid. getting spacing between k points
dk = 1/round(Nkpoints**(1/3))

# summarizing for each k point
sum_over_cv = np.zeros((Nkpoints), dtype = float)

for ik in range(Nkpoints):
    sum_over_cv[ik] = np.sum(abs(Akcv[ik, :, :]))
    
# scale it to make maximum value = 1
sum_over_cv = sum_over_cv / np.max(sum_over_cv)
    
    

def sphere_Akcv(radius, center):
    # Define sphere parameters
    # r = 2
    # center = (1, 2, 3)

    # Create sphere mesh
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    # in the np.meshgrid, the j variable is used to specify the number of points to use in the grid, so it is not related to the imaginary unit.
    x = center[0] + radius * np.sin(v) * np.cos(u)
    y = center[1] + radius * np.sin(v) * np.sin(u)
    z = center[2] + radius * np.cos(v)
    
    return x, y, z

# Create a 3D figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# plotting spheeres
for ik in range(Nkpoints):
    radius = sum_over_cv[ik] * dk 
    center = Kpoints_bse[ik]
    
    if put_everything_1BZ == True:
        for idir in range(3):
            if center[idir] < 0:
                center[idir] = center[idir] + 1
    # Plot the sphere
    
    if radius >= 1e-2:
        print(center, radius)
        x, y, z = sphere_Akcv(radius, center)
        ax.plot_surface(x, y, z, alpha=1,  color='red')
        
        # if plot_periodic_images == Tru e:
            
            # new_center = center + np.array([1,0,0])
            # x, y, z = sphere_Akcv(radius, new_center)
            # ax.plot_surface(x, y, z, alpha=1,  color='red')
            
            # new_center = center + np.array([0,1,0])
            # x, y, z = sphere_Akcv(radius, new_center)
            # ax.plot_surface(x, y, z, alpha=1,  color='red')

            # new_center = center + np.array([0,0,1])
            # x, y, z = sphere_Akcv(radius, new_center)
            # ax.plot_surface(x, y, z, alpha=1,  color='red')
            

            # new_center = center + np.array([0,1,1])
            # x, y, z = sphere_Akcv(radius, new_center)
            # ax.plot_surface(x, y, z, alpha=1,  color='red')
            
            # new_center = center + np.array([1,0,1])
            # x, y, z = sphere_Akcv(radius, new_center)
            # ax.plot_surface(x, y, z, alpha=1,  color='red')
            
            # new_center = center + np.array([1,1,0])
            # x, y, z = sphere_Akcv(radius, new_center)
            # ax.plot_surface(x, y, z, alpha=1,  color='red')
            
            # new_center = center + np.array([1,1,1])
            # x, y, z = sphere_Akcv(radius, new_center)
            # ax.plot_surface(x, y, z, alpha=1,  color='red')
            
        
# write letters
ax.text(0.5, 0.5, 0.5, "R", color="blue", fontsize=20)
ax.text(0.5, 0.5, 0.0, "M", color="blue", fontsize=20)
ax.text(0.5, 0.0, 0.0, "X", color="blue", fontsize=20)
ax.text(0.0, 0.0, 0.0, r"$\Gamma$", color="blue", fontsize=20)
        
        
# # Define points for line
# x, y, z = np.array([0, 1]), np.array([0, 0]), np.array([0, 0])
# ax.plot(x, y, z, color='black')

# x, y, z = np.array([0, 0]), np.array([0, 1]), np.array([0, 0])
# ax.plot(x, y, z, color='black')

# Define the corners of the cube
corners = np.array([[0, 0, 0],
                    [0, 0, 1.0],
                    [0, 1, 0],
                    [1, 0, 0],
                    [1, 1, 0],
                    [1, 0, 1],
                    [0, 1, 1],
                    [1, 1, 1]])

corners -= np.array([1/2, 1/2, 1/2])

# Define the edges of the cube
edges = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4), (2, 6), (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)]
edges = [(0, 1), (0, 2), (0, 3), (1, 5), (1,6), (2, 4), (2,6), (3, 4), (3,5), (4, 7), (5, 7), (6,7)]

# 
for edge in edges:
    corner_from = corners[edge[0]]
    corner_to = corners[edge[1]]
    x, y, z = np.array([corner_from[0], corner_to[0]]), np.array([corner_from[1], corner_to[1]]), np.array([corner_from[2], corner_to[2]])
    ax.plot(x, y, z, color='black')
    
# # Draw the cube edges
# edge_collection = Line3DCollection([corners[edge] for edge in edges], colors='black')
# ax.add_collection(edge_collection)



# Set axis labels
ax.set_xlabel('kx')
ax.set_ylabel('ky')
ax.set_zlabel('kz')
lim = 0.75
ax.set_xlim3d(-lim, lim)
ax.set_ylim3d(-lim, lim)
ax.set_zlim3d(-lim, lim)

# Show the plot
plt.show()