import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shutil
import csv

def parse_xfoil_polar(filepath):
    """reads xfoil_polar data"""
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
                    'CDp' : float(cols[3]),
                    'CM' : float(cols[4]),
                    'Top_Xtr' : float(cols[5]),
                    'Bot_Xtr': float(cols[6])
                    
                    
                }
                data.append(entry)
            except ValueError:
                continue
    return data


def aggregate_results(input_folder, assignment_dir, param_variation):
    final_results_path = os.path.join(assignment_dir, f"{param_variation}_full_results_df.csv")
    all_data = []


    #name is like f"{airfoil_name}_polars_raw_{r_name}_T_{thick: .3f}_{U_inf_mph}_{RPM}.txt"
    # Cerchiamo tutti i file .txt nella cartella
    for filename in os.listdir(input_folder):
        if ("polars") in filename and filename.endswith(".txt"):
            varying_param = filename.split('_')[-1].split('.txt')[0]
            
            full_path = os.path.join(input_folder, filename)
            
            # Usiamo il tuo parser originale
            points = parse_xfoil_polar(full_path)
            
            for p in points:
                p[param_variation] = varying_param
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


def run_xfoil_Cl( filename, assignment_dir, Re,  xtr_top_val = 1,  Cl = 0.4, alpha = 5):
    # Nome del file temporaneo per i risultati


    xfoil_path = r"C:\copiarefiles\XFOIL\xfoil_windows\working_version\xfoil.exe"
    worker_path, xfoilname = os.path.split(xfoil_path)

    res_path =  os.path.join(assignment_dir, filename)

    temp_local_path = os.path.join(worker_path, filename)
    
    polars_acc_path = os.path.join(assignment_dir, 'temp_polars')
    final_dest_path = os.path.join(polars_acc_path, filename)
    
    if os.path.exists(temp_local_path): os.remove(temp_local_path)
    if os.path.exists(final_dest_path): os.remove(final_dest_path)

    if os.path.exists(res_path): os.remove(res_path)

        # Comandi da inviare a XFOIL
    commands = (
        "PLOP\n"
        "G\n"
        "\n"
        "naca 2310\n"
        "OPER\n"
        "MACH 0\n"
        f"VISC {int(Re)}\n"
        "ITER 200\n"
        "VPAR\n"
        "N 9\n"
        "xtr\n"
        f"{xtr_top_val}\n"
        "\n"
        "\n"
        
        "PACC\n"
        f"{filename}\n"
        "\n"
        f"CL {Cl}\n"
        "\n"
        "QUIT\n"
    )
   
    
    
    
    
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
                                    cwd=worker_path)
            
            print("XFOIL returncode:", proc.returncode)
            print("XFOIL stdout:\n", proc.stdout)
            print("XFOIL stderr:\n", proc.stderr)
           
            
            # 3. SPOSTAMENTO FILE
            if os.path.exists(temp_local_path):
                shutil.move(temp_local_path, final_dest_path)
                print(f"File {filename} generato e spostato in {polars_acc_path}")
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


def run_xfoil_alfa( filename, assignment_dir, Re,  xtr_top_val = 1,  Cl = 0.4, alpha = 0):
    # Nome del file temporaneo per i risultati


    xfoil_path = r"C:\copiarefiles\XFOIL\xfoil_windows\working_version\xfoil.exe"
    worker_path, xfoilname = os.path.split(xfoil_path)

    res_path =  os.path.join(assignment_dir, filename)

    temp_local_path = os.path.join(worker_path, filename)
    
    polars_acc_path = os.path.join(assignment_dir, 'temp_polars_Re')
    final_dest_path = os.path.join(polars_acc_path, filename)
    
    if os.path.exists(temp_local_path): os.remove(temp_local_path)
    if os.path.exists(final_dest_path): os.remove(final_dest_path)

    if os.path.exists(res_path): os.remove(res_path)

        # Comandi da inviare a XFOIL
    commands = (
    "PLOP\n"
    "G\n"
    "\n"
    "naca 2310\n"
    "OPER\n"
    "MACH 0\n"
    f"VISC {int(Re)}\n"
    "VPAR\n"
    "N 12\n"
    "xtr\n"
    f"{xtr_top_val}\n"
    "\n"
    "\n"
    # qui fai il calcolo
    f"ALFA {alpha}\n"
    
    # ORA che la soluzione è stata calcolata, salva
    "CPWR\n"
    f"{filename.replace('.txt','')}_Cp.txt\n"
    "DUMP\n"
    f"{filename.replace('.txt','')}_Dump.txt\n"
    "\n"
    "\n"
    "QUIT\n"
)
    
    
    
    
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
                                    cwd=worker_path)
            
            print("XFOIL returncode:", proc.returncode)
            print("XFOIL stdout:\n", proc.stdout)
            print("XFOIL stderr:\n", proc.stderr)
           
            
            # 3. SPOSTAMENTO FILE
            if os.path.exists(temp_local_path):
                shutil.move(temp_local_path, final_dest_path)
                print(f"File {filename} generato e spostato in {polars_acc_path}")
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


def top_transition_variation(assignment_dir, xtr_lst, Re = 800000, Cl = 0.4, alpha = None ):
    
    file_acc_path =  os.path.join(assignment_dir, 'temp_polars')
    if not alpha:
        alpha = 'no'
    if not Cl:
        Cl = 'no'
    for element in xtr_lst:
        element = round(element, 2)
        filename = f'polars_2310_Re{Re}_Cl{Cl}_alpha{alpha}_xtrtop{element}.txt'
        run_xfoil_Cl(filename, assignment_dir, element)
    aggregate_results(file_acc_path, assignment_dir, 'xtr_top')


def laminar_separation_bubble_resarch(assignment_dir, Re_lst, xtr =1,  Cl = None, alpha = 0 ):
    
    file_acc_path =  os.path.join(assignment_dir, 'temp_polars_Re')
    if not alpha:
        alpha = 'no'
    if not Cl:
        Cl = 'no'
    for element in Re_lst:
        element = round(element)
        filename = f'polars_2310_xtr{xtr}_Cl{Cl}_alpha{alpha}_Re{element}.txt'
        run_xfoil_alfa(filename, assignment_dir, element)
    aggregate_results(file_acc_path, assignment_dir, 'Reynolds')

def laminar_sep_bubble_plot():
    sep_bubble_path = r"C:\copiarefiles\aircraft_ady\assignments_material\assignment_2\temp_polars_Re"
    sep_bubble_df = aggregate_results(sep_bubble_path, assignment_dir, 'Reynolds')
    Dump_bubble = []

    Cp_bubble = []
    for root, dir, filetouples in os.walk(sep_bubble_path):   
        for file in filetouples:
        
            if '_Dump' in file:
                sub_df = pd.read_csv(os.path.join(root, file), delimiter= '\s+', comment='#', 
                                     names = ["s", "x", "y", "Ue_Vinf", "Dstar", "Theta", "Cf", "H"])
                sub_df['Reynolds'] = file.split('_')[-2]
                if len(sub_df)>1:
                    Dump_bubble.append(sub_df)
            elif '_Cp' in file:
                sub_df_cp = pd.read_csv(os.path.join(root, file), delimiter= '\s+', comment='#', skiprows=2,
                                     names = ["x", "y", "Cp"])
                sub_df_cp['Reynolds'] = file.split('_')[-2]
                if len(sub_df_cp)>1:
                    Cp_bubble.append(sub_df_cp)
                
            
            
        
        # concatenazione corretta:
    dump_bubble_df = pd.concat(Dump_bubble, ignore_index=True)
    cp_bubble_df = pd.concat(Cp_bubble, ignore_index=True)

    dump_bubble_df.to_csv(r"C:\copiarefiles\aircraft_ady\assignments_material\assignment_2\dump_files_sep_bubble.csv")
    cp_bubble_df.to_csv(r"C:\copiarefiles\aircraft_ady\assignments_material\assignment_2\cp_files_sep_bubble.csv")


    for Re_val ,group in cp_bubble_df.groupby('Reynolds'):
        #unfortunately zip() does not work, need to create a mask n  the other df
        dump_grouped = dump_bubble_df[dump_bubble_df['Reynolds'] == Re_val]
        if dump_grouped.empty:
            continue #IF NOT THE GRAPHS WILL SHIFT 

        #I saw that the coordinates are "suction" and pressure based on x increasing or not. this is my geometric sep
        i_le = group['x'].idxmin()  # leading edge, finds that (where touches 0) and use to divide

        suction_raw = group.loc[:i_le]   #remember, this is still a df. can use the same principles
        pressure_raw = group.loc[i_le:]
        suction = suction_raw.sort_values('x')
        pressure = pressure_raw.sort_values('x')

        #repeat for dump
        i_le_dump = dump_grouped['x'].idxmin()
        dump_upper_raw = dump_grouped.loc[:i_le_dump]
        dump_upper = dump_upper_raw.sort_values('x')
    
        fig, axs = plt.subplots(1, 2, figsize=(12, 4), sharex=True)

        # ================= Cp plot =================
        axs[0].plot(
            suction['x'],
            suction['Cp'],
            '-o',  color = 'cyan',
            markersize=2,
            label='Cp suction side'
        )
        axs[0].plot(
            pressure['x'],
            pressure['Cp'], 
            '-o', color = 'orange',
            markersize=2,
            label='Cp pressure side'
        )

        axs[0].invert_yaxis()
        axs[0].set_xlabel('x/c')
        axs[0].set_ylabel('Cp')
        axs[0].set_title(f'Cp distribution at {Re_val})')
        axs[0].legend()
        axs[0].grid(True)

        # ================= Cf plot =================
        axs[1].plot(
            dump_upper['x'],
            dump_upper['Cf'],
            '-o',
            markersize=2,
            label='Cf upper'
        )

        # Highlight separated region (bubble)
        bubble = dump_upper['Cf'] <= 0
        axs[1].scatter(
            dump_upper.loc[bubble, 'x'],
            dump_upper.loc[bubble, 'Cf'],
            color='red',
            s=10,
            label='Separated'
        )

        axs[1].axhline(0.0, color='k', linestyle='--', linewidth=0.8)
        axs[1].set_xlabel('x/c')
        axs[1].set_ylabel('Cf')
        axs[1].set_title('Skin friction Cf (upper side)')
        axs[1].grid(True)
        axs[1].legend()

        plt.savefig(os.path.join(r"C:\copiarefiles\aircraft_ady\assignments_material\assignment_2\figures", f'{Re_val}_bubble.png'))
        


def plot_x_foil_polars(input_path, plot = True):
    output_path = r"C:\copiarefiles\aircraft_ady\assignments_material\assignment_2\figures"
    df_raw = pd.read_csv(input_path,  skiprows = 10, sep = '\\s+')
    mask = df_raw['alpha'].str.contains('---')
    df = df_raw[~mask]
    df = df.apply(pd.to_numeric, errors = 'coerce')
    Cl = df['CL']
    Cd = df['CD']
    alpha = df['alpha']
    xtr_top = df['Top_Xtr']
    xtr_bot = df['Bot_Xtr']
    if plot:
        plt.plot(alpha, Cl, color = 'blue', label = 'Cl', marker='o', markersize=2)
        plt.plot(alpha, Cd, color = 'red', label = 'Cd', marker='o', markersize=2)
        plt.legend()
        plt.title(r'Polar: Cl - $\alpha$ vs Cd - $\alpha$')
        plt.xlabel(r'angle of attack $\alpha [\circ]$ ')
        plt.ylabel('Cl , Cd')
        plt.savefig(os.path.join(output_path, 'ClvsCd_Re800000_alpha_sm2e8pp5.png'))
        plt.close()
    

        plt.plot(alpha, xtr_bot, color = 'green', label = 'Xtr_bot', marker='o', markersize=2)
        plt.plot(alpha, xtr_top, color = 'lime', label = 'Xtr_top', marker='o', markersize=2)
        plt.legend()
        plt.title(r'Transition points with different - $\alpha$ ')
        plt.xlabel(r'angle of attack $\alpha [\circ]$ ')
        plt.ylabel('Xtr_top, Xtr_bot')
        plt.savefig(os.path.join(output_path, 'xtr_var_Re800000_alpha_sm2e8pp5.png'))
        plt.close()
    return Cl, Cd, alpha, xtr_top, xtr_bot



if __name__ == '__main__':
    Re_lst = np.linspace(100000, 700000, 100)
    assignment_dir = r"C:\copiarefiles\aircraft_ady\assignments_material\assignment_2"
   
    laminar_sep_bubble_plot()
    plot_x_foil_polars(r"C:\copiarefiles\aircraft_ady\assignments_material\assignment_2\s_2310_M0_Re800000_n9")

