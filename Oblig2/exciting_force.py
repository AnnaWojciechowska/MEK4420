import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.special
import warnings
#warnings.simplefilter("error")
#warnings.filterwarnings('error', category=ComplexWarning)

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

def exciting_force_comp(file_path, L, D):
    x_m, y_m, x_p, y_p,  nn_total_segments = parse_matlab_output(file_path)
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
     
    gg = np.zeros([nn_total_segments,nn_total_segments], dtype=complex)
    ss = np.zeros([nn_total_segments,nn_total_segments], dtype=complex)

    wave_number_array = np.arange(0.01, 2, 0.1)
    XX2_sum = np.empty(len(wave_number_array), dtype=complex)
    Haskind1_sum = np.empty(len(wave_number_array), dtype=complex)
    Haskind2_sum = np.empty(len(wave_number_array), dtype=complex)
    Froude_Krylov_sum = np.empty(len(wave_number_array), dtype=complex)
    
    for k_index, nu in enumerate(wave_number_array):
        phi0 = np.exp(nu*(by- complex(0, 1) * bx))
    
        for i in range(0,nn_total_segments):
            for j in range(0, nn_total_segments):
                xa = bx[j] - bx[i]
                yb = by[j] + by[i]
                zz = nu*(yb-complex(0,1)*xa)
                # expint exponential integral function
                f1 = -2*np.exp(zz)*(scipy.special.exp1(zz)+np.log(zz)-np.log(-zz))
                f2 = 2*np.pi*np.exp(zz)
                g2 = np.real(f1)+complex(0,1)*np.real(f2)

                #exciting force:
                xa1 = xg1[j] - bx[i]
                xa2 = xg2[j] - bx[i]
                ya1 = yg1[j] - by[i]
                ya2 = yg2[j] - by[i]
                ra1 = np.sqrt(xa1*xa1+ya1*ya1)
                ra2 = np.sqrt(xa2*xa2+ya2*ya2)
                g0 = (np.log(ra1) + np.log(ra2))*0.5
                #% all other terms with midpoint rule
                rb = np.sqrt(xa*xa+yb*yb)
                g1 = -np.log(rb)
                gg[i,j] = (g0+g1+g2)*ds[j]
                    
                #lhs
                arg0 = np.imag(np.log((x_m[j]-bx[i]+complex(0,1)*(y_m[j]- by[i]))/(x_p[j]-bx[i]+complex(0,1)*(y_p[j]-by[i]))))
                if j-i == 0:
                    arg0=-np.pi
                arg1 = np.imag(np.log((x_m[j]-bx[i]+complex(0,1)*(y_m[j]+by[i]))/(x_p[j]-bx[i]+complex(0,1)*(y_p[j]+by[i]))))
                help1 = (n1[j]* (np.imag(f1)+complex(0,1)*np.imag(f2)) +n2[j]* (np.real(f1)+complex(0,1)*np.real(f2)))*nu*ds[j]
                ss[i,j]=(arg0+arg1+help1)
                

        # computation for Haskind relations:
        rhs = np.dot(gg,n2)
        phi_2=np.linalg.solve(ss,rhs)   
        #Haskind1
        XX2_Haskind1 =(phi0 *n2 - phi_2*nu*(n2- n1*complex(0,1))*phi0)*ds
        Haskind1_sum[k_index] = np.abs(XX2_Haskind1.sum()*complex(0,1))
        
        #Haskind2
        #far amplitude in minus infinity
        FarAmplitude2 = (phi_2*(nu*n2 - nu*complex(0,1)*n1)-n2) * phi0*ds
        Haskind2_sum[k_index] = np.abs(FarAmplitude2.sum()*complex(0,1))

        # Froude Krylov
        Froude_Krylov_sum 
        # D is 1 for all of the goemetries
        F_K = L * np.exp(-nu*D) *np.sin(nu*L/2)/(nu*L/2)
        Froude_Krylov_sum[k_index]= np.abs(F_K.sum())
        
        
        # diffraction
        rhsD = - 2. * np.pi * phi0
        phi_diffraction = np.linalg.solve(ss, rhsD)
        #exciting force from pressure:
        XX2 = phi_diffraction*n2*ds
        XX2_sum[k_index] = np.abs(XX2.sum()*complex(0,1))
    
    return XX2_sum, Haskind1_sum, Haskind2_sum, Froude_Krylov_sum


data_dir_path = '/home/anna/annaCode/UiO/MEK4420/Oblig2_work/'
file_list = ['box10_1.dat','box2_1.dat','box1_1.dat', 'box1_10.dat']
#file_list = ['box2_1.dat']
aspect_ratio_strings = ['L:D 10:1', 'L:D 2:1', 'L:D 1:1', 'L:D 1:10']
aspect_ratios = [ (10,1), (2,1),(1,1), (1,10)]

fig, axs = plt.subplots(2, 2, figsize=(12, 8), tight_layout=True)

for i in range(0,4):
    X_pressure, Haskind1, Haskind2, Froude_Krylov =  exciting_force_comp(os.path.join(data_dir_path, file_list[i]), aspect_ratios[i][0], aspect_ratios[i][1])
    wave_number_array = np.arange(0.01, 2, 0.1)
   
    #axs[i//2, i%2].plot(wave_number_array, X_pressure, color = 'red', alpha = 0.5, label= "X2 from pressure")
    axs[i//2, i%2].plot(wave_number_array, X_pressure, label = "pressure integration")
    axs[i//2, i%2].plot(wave_number_array, Haskind1, '--', label = "Haskind1")
    axs[i//2, i%2].plot(wave_number_array, Haskind2,'*', label = "Haskind2")
    axs[i//2, i%2].plot(wave_number_array, Froude_Krylov, label = "Froude_Krylov")
    axs[i//2, i%2].set_title('comparision for exciting force for ' + aspect_ratio_strings[i])
 
    axs[i//2, i%2].set_title('comparision of computing exciting force for ' + aspect_ratio_strings[i])
    axs[i//2, i%2].set_ylabel(r'$\frac{\left|X_2\right|}{\rho g D}$',  rotation=0)
    axs[i//2, i%2].set_xlabel(r'$\frac{\omega^{2}D}{g}$')
    axs[i//2, i%2].grid(True)
    axs[i//2, i%2].legend()

plt.legend()
plt.show()
