#this file should be done to create a BEM solver

import os 
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import csv
import time
import shutil
from collections import defaultdict
from scipy.optimize import brentq
# BEM : union between momentum theory (with also induction factor in radial direction) and loads from blade elements



def parse_xfoil_polar(filepath):
    """Legge TUTTE le righe del file polare prodotto da XFOIL"""
    data = []
    if not os.path.exists(filepath):
        return data
        
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
    # In theory the xfoil out data are all the same, they have ---- and then data
    start_index = 0
    for i, line in enumerate(lines):
        if "------" in line:
            start_index = i + 1
            break
            
    # this happens if it doesn't find the header
    if start_index == 0 or start_index >= len(lines):
        print('no header found')
        return data

    for line in lines[start_index:]:
        cols = line.split()
        # make sure there are enough columns (if not it is not an actual col)
        if len(cols) >= 3: 
            try:
                entry = {
                    'alpha': float(cols[0]),
                    'CL': float(cols[1]),
                    'CD': float(cols[2]),
                    
                }
                data.append(entry)
            except ValueError:
                continue
    return data


def aggregate_results(input_folder, working_dir, airfoil_name):
    final_results_path = os.path.join(working_dir, f"{airfoil_name}_full_results_df.csv")
    all_data = []


    #name is like f"{airfoil_name}_polars_raw_{r_name}_T_{thick: .3f}_{U_inf_mph}_{RPM}.txt"
    # Cerchiamo tutti i file .txt nella cartella
    for filename in os.listdir(input_folder):
        if ("polars_raw") in filename and airfoil_name in filename and filename.endswith(".txt"):
            rpm = filename.split("_")[-1].split('.txt')[0] #very annoying, it contains the txt info 
            thick_ratio = filename.split("_")[-3]
            velocity = filename.split("_")[-2]
            station_idx = filename.split("_")[-5]
            full_path = os.path.join(input_folder, filename)
            
            # Usiamo il tuo parser originale
            points = parse_xfoil_polar(full_path)
            
            for p in points:
                p['station_idx'] = station_idx
                p['U_inf'] = velocity
                p['rpm'] = rpm
                p['thick_ratio'] = thick_ratio
                all_data.append(p)

    if all_data:
        keys = all_data[0].keys()
        with open(final_results_path, 'w', newline='') as f:
            dict_writer = csv.DictWriter(f, fieldnames=keys)
            dict_writer.writeheader()
            dict_writer.writerows(all_data)
        print(f"\nSuccesso! Creato file unico con {len(all_data)} punti.")
    else:
        print("\nNessun dato trovato da aggregare.")
    return all_data


def run_xfoil(airfoil_path, alphas_xfoil, Re, thick, filename, working_dir):
    # Nome del file temporaneo per i risultati

    # first of all, I want to delete all the files that were generated in the past run
    alpha_start, alpha_end, alpha_step = alphas_xfoil

    scale_factor = thick / 0.1171   #to find new thickenss, use the one from clark
    airfoil_dir, airfoil_name = os.path.split(airfoil_path)  #if they are in same folder as xfoil it is easier
    res_path =  os.path.join(working_dir, filename)

    temp_local_path = os.path.join(airfoil_dir, filename)
    final_dest_path = os.path.join(working_dir, filename)
    
    if os.path.exists(temp_local_path): os.remove(temp_local_path)
    if os.path.exists(final_dest_path): os.remove(final_dest_path)

    print(airfoil_path)
    if os.path.exists(res_path): os.remove(res_path)
    if Re>150000:
        # Comandi da inviare a XFOIL
        commands = (
        "PLOP\n"
        "G\n"
        "\n"
        f"LOAD {airfoil_name}\n"
        "GDES"
        f"TFAC\n"
        f"{scale_factor}\n"
        f"{scale_factor}\n"
        "PPAR\n"
        "N 260\n"
        "\n"
        "\n"
        "PANE\n"
        "\n"
        "OPER\n"
        "MACH 0\n"
        f"VISC {int(Re)}\n"
        "ITER 200\n"
        "VPAR\n"
        "N 9\n"
        "\n"
        "PACC\n"
        f"{filename}\n"
        "\n"
        f"ASEQ {alpha_start} {alpha_end} {alpha_step}\n"
        "\n"
        "QUIT\n"
    )
    else:

        #found on github on calculations Cd0 = 0.09 for Re = 100 000, in my case even lower so Cd0 = 0.01 ? not too bad

        commands = (
        "PLOP\n"
        "G\n"
        "\n"
        f"LOAD {airfoil_name}\n"
        "GDES"
        f"TFAC\n"
        f"{scale_factor}\n"
        f"{scale_factor}\n"
        "PPAR\n"
        "N 260\n"
        "\n"
        "\n"
        "PANE\n"
        "\n"
        "OPER\n"
        "MACH 0\n"
        "ITER 200\n"
        "PACC\n"
        f"{filename}\n"
        "\n"
        f"ASEQ {alpha_start} {alpha_end} {alpha_step}\n"
        "\n"
        "QUIT\n"
    )
        print(f'{Re}')
    
    xfoil_path = r"C:\copiarefiles\XFOIL\xfoil_windows\working_version\xfoil.exe"
    
    
    with open("xfoil_debug.inp", "w") as f:
        f.write(commands)
    f.close()
    

    with open("xfoil_debug.inp", "r") as f_input:
        try:
            proc = subprocess.run([xfoil_path],
                                    stdin=f_input,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    text=True,
                                    cwd=airfoil_dir)
           
            
            # 3. SPOSTAMENTO FILE
            if os.path.exists(temp_local_path):
                shutil.move(temp_local_path, final_dest_path)
                print(f"File {filename} generato e spostato in {working_dir}")
                return True
            else:
                print(f"XFOIL non ha generato il file per Re={Re}. Convergenza fallita?")
                return False

        except subprocess.TimeoutExpired:
            proc.kill()
            print(f"Timeout XFOIL per {filename}")
            return False
        except Exception as e:
            print(f"Errore: {e}")
            return False


def retrieve_prop_data(csv_path):

    dataset = pd.read_csv(csv_path, sep=',', engine='python')
    
    for col in ['CHORD', 'THICKNESS', 'TWIST', 'STATION']:
        dataset[col] = pd.to_numeric(dataset[col], errors='coerce')

    dataset = dataset.dropna(subset = 'STATION')
    chord = dataset['CHORD'] 
    thick_ratio = dataset['THICKNESS'] 
    theta = dataset['TWIST']
    station = dataset['STATION'] 
    return station, chord, thick_ratio, theta


def create_polar_for_airfoil(airfoil_path, nu):
    '''this script is used to calculate the polar for an airfoil, to pass the correct Cl and Cd
    Ideally, we would also like to add the induction factor a, a' to the calculation.
    Of course, this would mean acting iteratively (it is a mess) and Re would not change for the induced velocity
    (it is ALWAYS) much smaller than the flow velocity. therefore we assume it to be a 2nd order effect and not consider it'''
    # Loop di accumulo
    final_results = []
    #get the geometrical values for my real propeller
    path_propeller_geo = r"C:\copiarefiles\aircraft_ady\assignments_material\assignment_5\prop_airfoil_char.csv"
    station_inches, chord_dist_inches, thick_ratio, _ = retrieve_prop_data(path_propeller_geo)
    inches_to_m = 2.54 / 100
    station = station_inches * inches_to_m
    chord_dist = chord_dist_inches * inches_to_m
    path_test_campaign = r"C:\copiarefiles\aircraft_ady\assignments_material\assignment_5\full_measurement_campaign.csv"
    measurement_campaign = pd.read_csv(path_test_campaign)
    U_inf_series = measurement_campaign['V']
    RPM_series = measurement_campaign['RPM']
    results_path = r"C:\copiarefiles\aircraft_ady\assignments_material\assignment_5"
    file_acc_path = os.path.join(results_path, 'temp_polars')

    if not os.path.exists(file_acc_path):
        os.makedirs(file_acc_path)

    dir, airfoil_name = os.path.split(airfoil_path)
    test_list = np.arange(0.02, 1.001, 0.02)
    alpha_lst = [-4, 12, 0.2]  # same pattern as AS in xfoil
    print(test_list)
    for i in range(len(station)):

        for U_inf_mph, RPM in zip(U_inf_series, RPM_series):
            
            U_inf = U_inf_mph * 0.447  #m/s
            U_inf = round(U_inf, 3)
            Omega = RPM / 60 * 2*np.pi  #rad/s
            r = station.iloc[i]  #this is the value already in meters, from the file
            #note there will be an error for the static data. since they don't use VTOT therec
            r_name = station_inches.iloc[i]
            thick = thick_ratio.iloc[i]
            filename = f"{airfoil_name}_polars_raw_{r_name}_T_{thick: .3f}_{U_inf_mph}_{RPM}.txt" 
            c = chord_dist.iloc[i]  #also in this case I have the chord, but I don't need it apart from the Re
            
            w = Omega  *  r
            v = np.sqrt(U_inf**2 + w**2)  #here lies the approximation, no induced vels
            Re  = (c * 0.75 ) * v / nu  # use half the chord length for Re (literature)
            
            successful_analysis = run_xfoil(airfoil_path, alpha_lst, Re, thick , filename, file_acc_path)
            if successful_analysis:
                print("OK")
            else:
                print('FAILED')

    aggregate_results(file_acc_path, results_path, airfoil_name)


def find_nearest_station(df, station, U_inf, rpm):
    # filtra prima per U_inf e rpm (o anche solo rpm se vuoi)
    df_sub = df[
        (df["rpm"] == int(rpm)) &
        (np.isclose(df["U_inf"], U_inf, atol=0.5))
    ]

    if df_sub.empty:
        raise ValueError(f"Nessuna polare per U_inf={U_inf}, rpm={rpm}")

    # distanza radiale
    idx = (df_sub["station_idx"] - station).abs().idxmin()
    return df_sub.loc[idx, "station_idx"]


def interpolate_Cl_Cd_curves_viterna( alpha_real, alpha_data, Cl_data, Cd_data, AR):
    
    Cd0_Clarky = 0.023
    #regression or corrections, remember that interp automatically manages the case where the data is already there 
    if alpha_real > np.min(alpha_data) and alpha_real < np.max(alpha_data) : #it is "in the curve"
        Cl_real = np.interp(alpha_real, alpha_data, Cl_data) #work with linear interp, since data is assumed to be "perfect" (from XFOIL)
        Cd_real = np.interp(alpha_real, alpha_data, Cd_data)
    else:
        #params for viterna
        Cd_max = 1.11 + 0.018* AR
        idx_stall = alpha_data.argmax()
        Cl_stall = Cl_data[idx_stall]
        Cd_stall = Cd_data[idx_stall]
        alpha_max = alpha_data[idx_stall]
        alpha_max_rad = np.deg2rad(alpha_max)

        alpha_real_rad = np.deg2rad(alpha_real)
        A_1 = Cd_max/ 2
        B_1 = Cd_max
        A_2 = ( Cl_stall - Cd_max*np.sin(alpha_max_rad)*np.cos(alpha_max_rad)) *np.sin(alpha_max_rad)/(np.cos(alpha_max_rad)**2)
        B_2 = (Cd_stall - Cd_max*(np.sin(alpha_max_rad)**2))/ np.cos(alpha_max_rad)
        #here there is also an approximation, in this case we assume that the highest alpha is the stall (not necessarily true)
        
        
        Cl_real = A_1 * np.sin(2 * alpha_real_rad) + A_2 * np.cos(alpha_real_rad)
        Cd_real = B_1 * np.sin(alpha_real_rad)**2 + B_2 *np.cos(alpha_real_rad)
    
    return Cl_real, Cd_real
    

def PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, NBlades, phi):
    eps = 1e6
    temp_tip = -(NBlades/2) * (tipradius_R-r_R)/(r_R * np.sin(phi) + eps)
    temp_tip = np.clip(temp_tip, -50, 50)
    Ftip = 2/np.pi * np.arccos(np.exp(temp_tip)) 

    temp_root = NBlades/2 * (rootradius_R-r_R)/(r_R * np.sin(phi) + eps)
    temp_root = np.clip(temp_root, -50, 50)
    Froot = 2/np.pi * np.arccos(np.exp(temp_root))

    F = np.clip(Ftip*Froot, 1e-3, 1.0)
    return F, Ftip, Froot

def find_omega_for_thrust(T_target, U_inf, Nb, chord_dist, theta_dist, dR, AR, cl_cd_df, Omega_min=1, Omega_max=3000):
    def f(Omega):
        # We only need the Thrust index (6th return value)
        results = solve_BEM(Nb, chord_dist, theta_dist, dR, AR, U_inf, Omega, cl_cd_df)
        thrust = results[5] 
        return thrust - T_target

    # Check if the target is within the bounds to avoid the Brentq crash
    f_min = f(Omega_min)
    f_max = f(Omega_max)

    if f_min * f_max > 0:
        # If both are negative, the prop is too small/slow to reach target thrust
        # If both are positive, it's already producing too much thrust
        print(f"Warning: Target thrust {T_target}N not reachable at U_inf={U_inf:.2f} m/s. Returning boundary.")
        return Omega_max if f_max < 0 else Omega_min

    Omega_solution = brentq(f, Omega_min, Omega_max)
    return Omega_solution


def solve_BEM(Nb, chord_dist, theta_dist, dR, AR, U_inf, Omega, cl_cd_df):
    #remember, the Cl works with alfa, but the alfa and phi are NOT the same angle
    U_inf_kph = max( U_inf / 0.44704, 0.01)
    iter = 0
    n = Omega / (2*np.pi)
    inches_to_m = 2.54 / 100
    a_prime = 0.01 * np.ones(n_segments) #azimuthal induction factor. assume there is only one value?
    a_prime_comp = np.zeros(n_segments)
    a = 0.1 * np.ones(n_segments)   #axial induction factor 
    a_comp = np.zeros(n_segments)
    
    hub_radius  = dR[0]
    prop_radius = dR[-1]
    rho, nu = get_air_quantities()
    prop_diam = prop_radius*2
    A_r = np.pi * (prop_diam/2) **2 #Area of the rotor (assuming all blade)

    rpm = int(round(Omega * 60 / (2*np.pi)))
    
    # T = 2 * rho *(U_inf - U_r)*U_r*A_r   # it is the same as doing the C_t later
    # C_t_2 = T / (0.5 * rho * A_r * U_inf**2)

    tol = 0.0001 #this breaks the loop
    #once the estimation for the C_t is done, we can calculate the velocities and the loads on the blade

    #important: instead of converging all together, the idea is that cross interactions
    #are not important, this means we can converge annuli separately


    #filter the df for the operating point, with the closest that we see
    available_rpms = cl_cd_df['rpm'].unique()
    closest_rpm = available_rpms[np.argmin(np.abs(available_rpms - rpm))]
    df_rpm = cl_cd_df[cl_cd_df['rpm'] == closest_rpm]
    available_uinf = df_rpm['U_inf'].unique()
    closest_uinf = available_uinf[np.argmin(np.abs(available_uinf - U_inf_kph))]
    df_op_point = df_rpm[df_rpm['U_inf'] == closest_uinf].copy()
    if df_op_point.empty:
        #means that maybe the we didn't have the analysis
        print(f'Warning: No analysis present for the case U_inf = {U_inf_kph}, rpm = {rpm}')
        return 0,0,0,0,0,0,0,0
    Thrust = 0
    Torque = 0

   
    converged = False
    Thrust = 0
    Torque = 0
    Power = 0
    for i in range(0, len(dR)):  #we want it to start at 1 since the calculation is in SEGMENTS
        iter = 0
        dr = dR[i] - dR[i-1] if i > 0 else dR[1] - dR[0]
        #remember to adjust the code so that 
        r = dR[i] #find the infinitesimal radious
        c = chord_dist[i]
        W_r = Omega * r* ( 1 - a_prime[i]) #find the velocity in the azimuthal direction
           
        U_r = U_inf * (1 + a[i])
        
        V_r = np.sqrt(W_r**2 + U_r**2 + 1e-9)  #this will be used to calculate the CD
        phi = np.atan2(U_r, W_r + 1e-19)

        #now we search for the polars (outside of the while, faster) based on the closest to station
        station_req = r / inches_to_m
        idx = (df_op_point['station_idx'] - station_req).abs().idxmin()
        station_use = df_op_point.loc[idx, "station_idx"]  #we select the correct index of the df
        df_polar = df_op_point[df_op_point['station_idx'] == station_use] #mask
        
        alpha_data = df_polar['alpha'].values
        Cl_data = df_polar['CL'].values
        Cd_data = df_polar['CD'].values

        #now we need to handle the inviscid case (Re was too low for calculation)
        if (Cd_data == 0).all():
            #use the simple version
            Cd0_Clarky = 0.023
            Cd_data = Cd0_Clarky + 0.045*(Cl_data**2)   #the classic formula
        
        #now we enter the loop 
        iter = 0
        converged = False
        while iter < n_iterations and converged == False:

            
            
            search_dict = {
                'station_idx' : r / inches_to_m,
                'rpm' : rpm,
                'U_inf_kph' : U_inf_kph
            }
    
            theta_rad = np.deg2rad(theta_dist[i])
            alfa_rad =  theta_rad - phi # if not would be a turbine, theta is pitch
            alfa_search = np.rad2deg(alfa_rad)
          

            Cl, Cd = interpolate_Cl_Cd_curves_viterna(alfa_search, alpha_data, Cl_data, Cd_data, AR)

               
            Lift = Cl * 0.5 * rho * V_r**2 * c  #per unit area? or not?
            Drag = Cd * 0.5 * rho * V_r**2 * c
            
            #just for the formula for a'
            

            #compute the new forces (the new Ct, Ctradial)
            dF_axial = Lift*np.cos(phi) - Drag*np.sin(phi)    #propeller ady (with plus) REMEMBER, THE FORCES ARE THE SUM OF BLADES
            dF_azimuthal = Lift*np.sin(phi) + Drag*np.cos(phi)
            dT  = Nb * dF_axial   #not only infinitesimal force, but also radius (needs to be multiplied)
            dQ = Nb*dF_azimuthal * r
            # Ct locale basato sulla sezione (annulus)
            
            F, _, _ = PrandtlTipRootCorrection(r/prop_radius, hub_radius/prop_radius, 1.0, Nb, phi)

            F = max(F, 1e-4)  # sicurezza numerica
            # now we have a difference, can use formulas from rotoe/wake or course. we will use course
            #the one from the course are from mass flow

            a_guess = dT / (2 * 2*np.pi * r * rho * U_inf**2 * (1 + a[i]) * F )

            a_prime_guess = dQ / (2*2 *np.pi * r**3 * Omega   * rho * U_inf *  (1+a[i])* F)

            
    

            a_new = 0.7*a[i] + 0.3* a_guess  #relax
            
        
            a_prime_new = 0.7*a_prime[i] + 0.3 * a_prime_guess

            a_new = np.clip(a_new, -0.2, 0.95)
            a_prime_new = np.clip(a_prime_new, -1.0, 1.0)
        
            a_comp[i] = a[i]
            a_prime_comp[i] = a_prime[i]
            a[i] = a_new
            a_prime[i] = a_prime_new
            iter +=1
            if abs(a_comp[i] - a[i]) < tol and abs(a_prime_comp[i] - a_prime[i])< tol:
                converged = True

        
        Thrust += dT * dr
        Torque += dQ * dr
        
        iter = 0
        
    Ct = Thrust / (rho * n**2 * prop_diam**4)
    Cq = Torque / (rho * n**2 * prop_diam**5)
    Cp = Cq * 2 *np.pi

    Power = Torque * Omega
    J = U_inf / (n*prop_diam)
    eta = (J * Ct) / (2 * np.pi * Cq)
    return Ct, Cq, Cp, eta, J, Thrust, Torque, Power


def prop_definition():

    #prop definition based on the 22x10 PROPELLER
    #basde on the definition, the first "airfoil" is here 3.1329, the first measured station (hub transition at 3.13)
    inches_to_m = 2.54 / 100
    n_segments = 15
    

    hub_rad = 3.13 * inches_to_m
    prop_r = 11 * inches_to_m  #CHANGE with propeller!
    blade_b = prop_r - hub_rad
    airfoil_path =  r"C:\copiarefiles\XFOIL\xfoil_windows\working_version\Clark_Y.dat"

    path_propeller_geo = r"C:\copiarefiles\aircraft_ady\assignments_material\assignment_5\prop_airfoil_char.csv"
    st_raw, ch_raw, _, th_raw = retrieve_prop_data(path_propeller_geo)

    st_raw_m = st_raw * inches_to_m
    ch_raw_m = ch_raw * inches_to_m

    dR = np.linspace(hub_rad, prop_r, n_segments)

    chord_interp = np.interp(dR, st_raw_m, ch_raw_m)
    theta_interp = np.interp(dR, st_raw_m, th_raw)
    min_chord = chord_interp.min()
    max_chord = chord_interp.max()
    AR = blade_b**2 / ((max_chord + min_chord) / 2 * blade_b)  #assuming a trapezoid (not so important) b^2 / S

    dir, airfoil_name_raw = os.path.split(airfoil_path)
    airfoil_name = airfoil_name_raw.split('.')[0]
    Cd0 = 0.01 
    Nb = 2 #number of blades
    
    
    
    return airfoil_name, Nb, prop_r, Cd0, chord_interp, theta_interp, dR, AR, n_segments 


def get_air_quantities():
    rho = 1.225
    nu = 1.4e-5
    return rho,nu


if __name__ == '__main__':

    airfoil_name, Nb, _, Cd0, chord_dist, theta_dist, dR, AR, n_segments = prop_definition()
    n_iterations = 50
    cl_cd_df_path = r"C:\copiarefiles\aircraft_ady\assignments_material\assignment_5\Clark_Y_full_results_df.csv"
    cl_cd_df = pd.read_csv(cl_cd_df_path)

    

    test_campaign_path = r"C:\copiarefiles\aircraft_ady\assignments_material\assignment_5\full_measurement_campaign.csv"

    full_campaign = pd.read_csv(test_campaign_path)
    
    M_inc = [0.09, 0.17, 0.18, 0.26, 0.27, 0.35, 0.44]  #last one should be incompressible limit

    incompressible_mask = np.zeros(len(full_campaign), dtype=bool)

    for M in M_inc:
        incompressible_mask |= np.isclose(
        full_campaign['Mach'].values,
        M,
        rtol=0,
        atol=1e-4
    )  #check only the ones that I have
        
    incompressible_campaign = full_campaign[incompressible_mask]  #this is the one that should work 
    compressible_campaign = full_campaign[~incompressible_mask] #this not

    for name, group in incompressible_campaign.groupby('Mach'):
        print(f'we are at Mach: {name}')
        J_ratio = group['J']
        Ct_bem = []
        Cp_bem = []
        eta_bem = []
        Cq_bem = []
        U_inf_array = group['V'].values
        
        perf = group['Pe'].values
        rpm_array = group['RPM'].values
        Ct_exp = group['Ct']
        Cp_exp = group['Cp']
        Cq_exp = Cp_exp / (2*np.pi)
        P_exp = group['PWR(W)']
        T_exp = group['Thrust(N)']

        eta_lst = []
        #scale_lst = [0.8, 0.9, 1.0, 1.1, 1.2]
        scale_lst = [1]  #use this to produce the result without varying the geometry
        if len(scale_lst) >1:
            for scale in scale_lst:
                
                J_curve = []
                eta_curve = []
                Ct_curve = []
                dR_new = dR * scale
                print(f'we are working with a radius of r = {dR_new}')
                for i in range(len(U_inf_array)):
                    if i%10 ==0:
                        print(f'working with velocity n = {i} of a total of {len(U_inf_array)}')
                    U_inf_raw = float(U_inf_array[i])
            
                    U_inf =max(U_inf_raw, 0.01)
                    rpm = rpm_array[i]
                    Omega = find_omega_for_thrust(0.08, U_inf, Nb, chord_dist, theta_dist, dR_new, AR, cl_cd_df)
                    Ct, Cq, Cp, eta, J, Thrust, Torque, Power = solve_BEM(Nb, chord_dist, theta_dist, dR_new, AR, U_inf, Omega, cl_cd_df )

                    eta_lst.append(eta)
                    J_curve.append(J)
                    eta_curve.append(eta)
                    Ct_curve.append(Ct)
                plt.plot(J_curve, eta_curve, marker='o', label=f'Scale {scale} (D={2*dR_new[-1] * 2}m)')


                # -----------------------------
                # LABEL & TITOLO
                # -----------------------------
            plt.xlabel('diameter', fontsize=12)
            plt.ylabel(' efficiency', fontsize=12)
            plt.title('Variation of efficiency with diameter ', fontsize=14, fontweight='bold')

                # -----------------------------
                # GRIGLIA
                # -----------------------------
            plt.grid(
                    which='both',
                    linestyle='--',
                    linewidth=0.7,
                    alpha=0.6
                )

                # -----------------------------
                # LEGENDA
                # -----------------------------
            plt.legend(
                    fontsize=11,
                    loc='best',
                    frameon=True,
                    edgecolor='black'
                )

                # -----------------------------
                # TICKS
                # -----------------------------
            plt.tick_params(
                    axis='both',
                    which='major',
                    labelsize=11,
                    direction='in'
                )

                # -----------------------------
                # MARGINI
                # -----------------------------
            plt.tight_layout()
            plt.show()

        else:
            # 1. Setup Geometry
            scale = scale_lst[0]
            dR_new = dR * scale
            # Important: Ensure D_scaled matches the experimental D
            prop_diam = dR_new[-1] * 2 
            
            J_bem, eta_bem, Ct_bem, Cp_bem , T_bem, P_bem = [], [], [], [], [], []
            
            # 2. Extract Experimental Vectors
            J_exp = group['J'].values
            Ct_exp = group['Ct'].values
            Cp_exp = group['Cp'].values
            eta_exp = group['Pe'].values
            
            print(f"Validating Mach {name} - {len(U_inf_array)} discrete points.")

            for i in range(len(U_inf_array)):
                # Use the EXACT same physical conditions as the test
                U_inf = float(U_inf_array[i])
                rpm_val = rpm_array[i]
                Omega_real = rpm_val * 2 * np.pi / 60
                
                # Run BEM for this specific point
                # Ensure Nb, chord_dist, etc., are the baseline ones
                Ct, Cq, Cp, eta, J, Thrust_bem, Torque_bem, Power_bem = solve_BEM(
                    Nb, chord_dist, theta_dist, dR_new, AR, U_inf, Omega_real, cl_cd_df
                )
                
                # If solve_BEM failed (returned 0s), skip it or append NaN
                if J == 0 and U_inf > 0:
                    J_bem.append(np.nan); eta_bem.append(np.nan)
                    Ct_bem.append(np.nan); Cp_bem.append(np.nan)
                    T_bem.append(np.nan); P_bem.append(np.nan)
                else:
                    J_bem.append(J)
                    eta_bem.append(eta)
                    Ct_bem.append(Ct)
                    Cp_bem.append(Cp)
                    T_bem.append(Torque_bem)
                    P_bem.append(Power_bem)

            # --- PLOTTING ---
            fig, axs = plt.subplots(1, 3, figsize=(18, 5))
            fig.suptitle(f'Point-to-Point Validation: Mach {name}', fontsize=14)

            metrics = [
                (T_bem, T_exp, '$C_t$', 'blue'),
                (P_bem, P_exp, '$C_p$', 'green'),
                (eta_bem, eta_exp, '$\eta$', 'black')
            ]

            for ax, (bem, exp, label, col) in zip(axs, metrics):
                # We use markers for both to see the X-axis alignment
                ax.plot(J_bem, bem, 'o-', label='BEM Result', color=col, alpha=0.7)
                ax.scatter(J_exp, exp, marker='x', color='red', label='Experimental', s=50, zorder=3)
                
                ax.set_xlabel('Advance Ratio $J$')
                ax.set_ylabel(label)
                # Set X-limit to match the experiments exactly
                ax.set_xlim(min(J_exp)*0.9, max(J_exp)*1.1)
                ax.legend()
                ax.grid(True, alpha=0.3)

            plt.tight_layout()
            plt.show()