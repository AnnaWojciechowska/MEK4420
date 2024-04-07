import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.special

plt.rcParams['text.usetex'] = True


def parse_matlab_output(file_path):
    num_lines = sum(1 for line in open(file_path))

    x_m = np.empty(num_lines)
    y_m = np.empty(num_lines)
    x_p = np.empty(num_lines)
    y_p = np.empty(num_lines)
    
    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            data = line.split()
            x_m[i] = float(data[0])  
            y_m[i] = float(data[1])  
            x_p[i] = float(data[2])  
            y_p[i] = float(data[3])  
    return x_m, y_m, x_p, y_p, num_lines


def calc_force_coefficents(file_path):
    wave_number_array = np.arange(0.01, 2, 0.1)
    ff_22_sum = np.empty(len(wave_number_array), dtype=complex)
    sum_AM2 = np.empty(len(wave_number_array), dtype=complex)
    sum_AP2 = np.empty(len(wave_number_array), dtype=complex)
  
    x_m, y_m, x_p, y_p,  nn_total_segments = parse_matlab_output(file_path)
    #sgement_inddex =np.linspace(1,nn_total_segments,nn_total_segments)
    
    dx_array=x_p - x_m
    dy_array=y_p - y_m   
    
    #.^2 element wise squaring in matlab
    ds = ((dx_array)**2+(dy_array)**2)**(1/2)
    # bx and by are midpoints
    bx = 0.5*(x_m + x_p)
    by = 0.5*(y_m + y_p)
    n1 = -(y_p - y_m)/ds
    n2 = (x_p - x_m)/ds
    
    
    #points for Gauss integration on each segment
    xg1 = -0.5 * dx_array/np.sqrt(3)+bx
    xg2 = 0.5  * dx_array/np.sqrt(3)+bx
    yg1 = -0.5 * dy_array/np.sqrt(3)+by
    yg2 = 0.5  * dy_array/np.sqrt(3)+by
    
    
    for k_index, nu in enumerate(wave_number_array):

        phi0 = np.exp(nu*(by- complex(0, 1) * bx))
        phi0n = nu * (n2 - complex(0,1) * n1)* phi0
    
        gg = np.zeros([nn_total_segments,nn_total_segments], dtype=complex)
        ss = np.zeros([nn_total_segments,nn_total_segments], dtype=complex)
        ss_radiation = np.zeros([nn_total_segments,nn_total_segments], dtype=complex)
    
    
        for i in range(0,nn_total_segments):
            for j in range(0, nn_total_segments):
            #rhs, log(r) term with 2pts Gauss quadrature
                xa1 = xg1[j] - bx[i]
                xa2 = xg2[j] - bx[i]
                ya1 = yg1[j] - by[i]
                ya2 = yg2[j] - by[i]
                ra1 = np.sqrt(xa1*xa1+ya1*ya1)
                ra2 = np.sqrt(xa2*xa2+ya2*ya2)
                g0 = (np.log(ra1) + np.log(ra2))*0.5
                #% all other terms with midpoint rule
                xa = bx[j] - bx[i]
                yb = by[j] + by[i]
                rb = np.sqrt(xa*xa+yb*yb)
                g1 = -np.log(rb)
                zz = nu*(yb-complex(0,1)*xa)
                # expint exponential integral function
                f1 = -2*np.exp(zz)*(scipy.special.exp1(zz)+np.log(zz)-np.log(-zz))
                f2 = 2*np.pi*np.exp(zz)
                g2 = np.real(f1)+complex(0,1)*np.real(f2)
                gg[i,j] = (g0+g1+g2)*ds[j]
                #lhs
                arg0 = np.imag(np.log((x_m[j]-bx[i]+complex(0,1)*(y_m[j]- by[i]))/(x_p[j]-bx[i]+complex(0,1)*(y_p[j]-by[i]))))
                if j-i == 0:
                    arg0 =np.pi
                arg1 = np.imag(np.log((x_m[j]-bx[i]+complex(0,1)*(y_m[j]+by[i]))/(x_p[j]-bx[i]+complex(0,1)*(y_p[j]+by[i]))))
                help1 = (n1[j]* (np.imag(f1)+complex(0,1)*np.imag(f2)) +n2[j]* (np.real(f1)+complex(0,1)*np.real(f2)))*nu*ds[j]
                ss[i,j]=(arg0+arg1+help1)
                # radiation problem
                arg0_radiation = np.imag(np.log((x_m[j]-bx[i]+complex(0,1)*(y_m[j]- by[i]))/(x_p[j]-bx[i]+complex(0,1)*(y_p[j]-by[i]))))
                if j-i == 0:
                    arg0_radiation = -np.pi
                arg1_radiation = np.imag(np.log((x_m[j]-bx[i]+complex(0,1)*(y_m[j]+by[i]))/(x_p[j]-bx[i]+complex(0,1)*(y_p[j]+by[i]))))
                help1_radiation = (n1[j]* (np.imag(f1)+complex(0,1)*np.imag(f2)) +n2[j]* (np.real(f1)+complex(0,1)*np.real(f2)))*nu*ds[j]
                ss_radiation[i,j] = (arg0_radiation + arg1_radiation + help1_radiation)        
        rhs = np.dot(gg,phi0n)
        rhs_radiation = np.dot(gg,n2)
        phi_radiation = np.linalg.solve(ss_radiation, rhs_radiation)

        #calculating of force coefficents:
        ff22 = (phi_radiation * n2 * ds)
        ff_22_sum[k_index] = np.sum(ff22)
        
        # calculations to obtain damping coeffient b22 from energy and far field amplitude:
        # this is not same phi0 as at the beginig of script
        complex(0,1)
        phi0_energy = np.exp(nu * (by -  complex(0,1) * bx))
        AM2 =  complex(0,1) * (phi_radiation * (nu * n2 - nu *  complex(0,1) * n1) - n2) * phi0_energy * ds
        AP2 =  complex(0,1) * (phi_radiation * (nu * n2 + nu *  complex(0,1) * n1) - n2) * np.conj(phi0_energy) * ds
        sum_AM2[k_index] = np.sum(AM2)
        sum_AP2[k_index] = np.sum(AP2)
    return ff_22_sum, sum_AM2, sum_AP2

file_list = ['box10_1.dat','box2_1.dat','box1_1.dat', 'box1_10.dat']
aspect_ratio_strings = ['L:D 10:1', 'L:D 2:1', 'L:D 1:1', 'L:D 1:10']

data_dir_path = '/home/anna/annaCode/UiO/MEK4420/Oblig2_work/'


swell_wave_len = 150
k_swell_wave = (2 *np.pi)/swell_wave_len

fig, axs = plt.subplots(2, 2, figsize=(12, 8), tight_layout=True)

for i, geometry_data in enumerate(file_list):
    ff_22_sum,  sum_AM2, sum_AP2 =  calc_force_coefficents(os.path.join(data_dir_path, geometry_data))
    wave_number_array = np.arange(0.01, 2, 0.1)
    b22_coef = -1 * np.imag(ff_22_sum)

    b22_energy = 0.5*(np.abs(sum_AM2)**2 + np.abs(sum_AP2)**2)
    #axs[i//2, i%2].plot(wave_number_array, b22_coef, label = r'$\frac{b_{22}}{\rho D^{2}}$')
    axs[i//2, i%2].fill_between(wave_number_array,b22_coef, b22_energy, color = 'red', alpha = 0.5, label= "difference")

    axs[i//2, i%2].set_title('comparision for damping coeff b22 obtained in 2 different ways for geometry ' + aspect_ratio_strings[i])
    axs[i//2, i%2].set_xlabel(r'$\frac{\omega^{2}D}{g}$')
    #  axs[i//2, i%2].axvline(k_swell_wave, color='navy', alpha = 0.5, linestyle='--')
   
    axs[i//2, i%2].grid(True)
    axs[i//2, i%2].legend()

plt.legend()
plt.show()
