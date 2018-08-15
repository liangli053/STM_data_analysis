#!/usr/bin/python
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from Savitzky_Golay_filter import savitzky_golay
from Process_Orig_Data import process_orig_data
from scipy.interpolate import interp1d

def get_1st_dev(input_data):
    ##-----calculate the 1st derivative, input data is a n x 2 matrix-----##
    tot_row = input_data.shape[0]
    out_matrix = np.zeros( (math.floor((tot_row-1)/2),2), dtype=float )
    aa = np.arange(1, tot_row-1,2)
    bb = np.arange(0, tot_row,2)
    out_matrix[:,0] = input_data[aa, 0]
    out_matrix[:,1] = np.diff(input_data[bb,1])/np.diff(input_data[bb,0])
    return out_matrix;

def integral(x,y):
    ##-----compute the integral of y over x-----#
    runningSum = 0
    for i in range(1,len(x)):
        width = x[i] -x[i-1]
        height = (y[i] + y[i-1])/2
        area = width*height 
        runningSum += area
    return runningSum

##-----find the original experimental data-----##
orig_dir = os.getcwd() + "/raw_data"
with open(orig_dir + "/file_names.txt" , 'r') as text:
  file_list = text.read().splitlines()
text.close()

with open(orig_dir + "/Z_vals.txt" , 'r') as text:
  Z = text.read().splitlines()
  Z = np.array(Z, float)
text.close()

## L_min and L_max are deternined in a separate run, these are the x-axis limits for interpolation
L_min = -12.0
L_max = 18.0
L_Range = np.arange(L_min, L_max, 0.01)
interpolated_plateau_all = np.empty((0, len(L_Range)), float)
os.chdir(orig_dir)

plt.figure()
for input_file in file_list:
    is_flat = 0
    curr_data = process_orig_data(input_file)  ## convert to plain text
    curr_data.get_scaling_factor()  ## get the scaling factor
    curr_data.extract_useful_data()  ## get the experimental data points at the end of the file
    curr_data.extract_plateau()  ## find the upper plateau
    curr_data.get_L()  ## calculate the scanned range
    curr_data.get_smooth_data_lowess() ## smooth data using Lowess method
    lowess_offset = curr_data.smooth_data_lowess
    if max(lowess_offset[:,1]) - min(lowess_offset[:,1]) > 0.14:
      x_offset = curr_data.smooth_data_lowess[ curr_data.smooth_data_lowess[:,1].argmin(), 0] ## x(L) coornidate of the valley
    else:
        valley_point = int(curr_data.smooth_data_lowess.shape[0]/2)
        x_offset = curr_data.smooth_data_lowess[valley_point, 0] ## find the mid-point value
        is_flat = 1

    lowess_offset = lowess_offset - [x_offset, 0] ## align the smooth data

##-----interpolate smoothed data, and use smooth data to do integration-----##
    if is_flat == 0:
      F = interp1d(lowess_offset[:,0], lowess_offset[:,1], kind ='cubic')
      bbb = F(L_Range)
      interpolated_plateau_all = np.append(interpolated_plateau_all, np.array([bbb]), axis=0)
    else:
      bbb = np.ones_like(L_Range) * lowess_offset[:,1].mean()
      interpolated_plateau_all = np.append(interpolated_plateau_all, np.array([bbb]), axis=0)

##-----show the messy-looking original data-----##
    #X_Y_all = curr_data.useful_data[:, 1:3]
    #X_Y_all = X_Y_all - np.ones_like(X_Y_all) * X_Y_all[0]
    #L_all = curr_data.scaling_factor * np.sqrt(np.square(X_Y_all[:,0]) + np.square(X_Y_all[:,1]) ) *3
    #plt.plot(L_all,curr_data.useful_data[:,3], '-')
    plt.plot(curr_data.useful_data[:,3], '-')

plt.show()
print interpolated_plateau_all.shape

fig, axes = plt.subplots(2,2)
for i in range(0, interpolated_plateau_all.shape[0]):
  interpolated_plateau_all[i,:] = interpolated_plateau_all[i,:] - interpolated_plateau_all[i,0]   ## adjust original plateaus

scale_fac = 3600.0/30414.0
kz = scale_fac * interpolated_plateau_all
for i in range(0, kz.shape[0]):
    axes[0,0].plot(L_Range, kz[i,:])
axes[0,0].set_title('kz')
axes[0,0].set_xlabel('X (angstrom)')
axes[0,0].set_ylabel('vertical stiffness (N/m)')

##-----integrate interpolated_plateau_all over Z to get Fz, Fz has one less row than fz-----##
Fz = np.zeros_like(kz)
for i in range(0, Fz.shape[0]):
    print i
    for j in range(0, Fz.shape[1]):
        Fz[i, j] = integral(Z[0:i+1], kz[0:i+1,j])

Fz = Fz * 100.0
for i in range(0, Fz.shape[0]):
    axes[0,1].plot(L_Range, Fz[i,:])
axes[0,1].set_title('Fz')
axes[0,1].set_xlabel('X (angstrom)')
axes[0,1].set_ylabel('vertical force (pN)')

##-----integrate Fz over Z to get Uz, Uz has one less row than Fz-----##
Uz = np.zeros_like(Fz)
for i in range(0, Uz.shape[0]):
    print i
    for j in range(0, Uz.shape[1]):
        Uz[i, j] = integral(Z[0:i+1], Fz[0:i+1,j] )

Uz = Uz/1.602
for i in range(0, Uz.shape[0]):
    axes[1,0].plot(L_Range, Uz[i,:])
axes[1,0].set_title('Uz')
axes[1,0].set_xlabel('X (angstrom)')
axes[1,0].set_ylabel('potential (meV)')

##-----derivative dUz/dL-----##
Uz = Uz * 1.602
for i in range(0, Uz.shape[0]):
    input_data = np.concatenate((L_Range[:, np.newaxis], Uz[i,:][:, np.newaxis]), axis = 1)
    aaa = get_1st_dev(input_data)
    axes[1,1].plot(aaa[:,0], aaa[:,1])
axes[1,1].set_title('Fx')
axes[1,1].set_xlabel('X (angstrom)')
axes[1,1].set_ylabel('lateral force (pN)')

mid_point = int(abs(L_min) / 0.01)
delta_f_2pts = interpolated_plateau_all[:, [0, mid_point]] ## record the delta_f values of each scan when L=left and L = mid
Fz_2pts = Fz[:, [0, mid_point]] ## record the Fz values of each scan when L=left and L = mid
Uz_2pts = Uz[:, [0, mid_point]] ## record the Uz values of each scan when L=left and L = mid
delta_f_Fz_Uz = np.concatenate((delta_f_2pts, Fz_2pts, Uz_2pts), axis = 1)
write_to = open("delta_f_Fz_Uz.txt", 'w')
for i in range(0, delta_f_Fz_Uz.shape[0]):
   write_to.write("  ".join(str("%8.3f" %e) for e in delta_f_Fz_Uz[i,:]) + '\n')
write_to.close()

plt.show()
