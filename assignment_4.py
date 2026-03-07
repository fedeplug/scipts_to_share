import os 
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shutil
import csv
import re

class Wing():
    #the wing is the most simple component in this case. not a real parametrization but the one given in assignment
    def __init__(self, b_1, b_2, root_c, tip_c, dihedral):
        self.b_1 = b_1
        self.b_2 = b_2 - b_1
        self.root_c = root_c
        self.kink_c = root_c
        self.tip_c = tip_c
        self.dihedral = dihedral
    
    def MAC(self):
        S_rect = self.b_1 * self.root_c
        S_trap = self.b_2 * (self.kink_c + self.tip_c) / 2
        mac_2 = 2/3 * (self.kink_c + self.tip_c - ((self.kink_c*self.tip_c)/(self.kink_c + self.tip_c)))
        return ((S_rect * self.root_c) + S_trap*mac_2) / (S_rect + S_trap)
    
    def calculate_sweep_deg(self):
        #in this case the sweep is only the angle between the final two points 
        sweep_rad = np.atan2((self.root_c - self.tip_c) ,self.b_2)
        return np.rad2deg(sweep_rad)
        
    def ady_centre_coords(self):
        S_ref =  self.b_1 * self.root_c + self.b_2 * (self.kink_c + self.tip_c) / 2
        mac = self.MAC()
        x_ref = 0.25 * mac
        y_ref = 0
        z_ref = 0
        return S_ref, x_ref, y_ref, z_ref
    


class Winglet():
    #the original winglet analyzed is taken from Airbus A320 (modeled in thesis : https://webthesis.biblio.polito.it/14640/1/tesi.pdf)
    #tha A320 in the end has the double "fence" config but they analyze what would change, let's use this normalized
    #original b = 2106.3 , root_c = 1600 , root_c_w = 1590 but use 1600 as per assignment
    #  tip_c = 435.0 (so we have the taper), delta Le = 67.1, delta Te = 134.1 so basically a bit of taper I think
    #cant angle not specified but I will take a straight one to start, so 90 deg, as alludes in the thesis
    #the twist is the only angle not specified. in paper I found -0.6 as I.C changed until 20 in study.
    # since I also don't have a 3d geom, I will take it as "added camber" so i will use -2
    def __init__(self, b, root_c, tip_c, delta_le, cant_angle, tip_twist ):
        self.b = b
        self.root_c = root_c
        self.taper = tip_c / root_c  #since I will pass the original "A32 normalized" data
        self.delta_le = delta_le  #IMPORTANT. to achieve this , delta_le will be measured angle at LE (also counting Taper). 
        #so we can calculate Te from the chord at tip that is fixed 
        self.cant_angle = cant_angle
        self.root_twist = 0 #since my plane will always start from 0 config (also in the thesis)
        self.tip_twist = tip_twist

        self.tip_c = root_c * self.taper  #circular importing but don't care, the issue is in A320 data
        self.tip_1 = 435 / 1600
        
    
    def build_geom(self, n_segments, start_point, wing_sweep_LE_deg):  # I define a lot of segments to avoid the 
        #90 deg sharp turn (should not make a difference I believe (no interference turbulence, AVL) but slightly less Cl, less surface exp)
        transition_length =  0.2  #random val we will see
        winglet_file_path = r"C:\copiarefiles\aircraft_ady\ADY_AVL\winglet.txt"
    

        with open(winglet_file_path, 'w') as file_handle:
            file_handle.write(f"#-------------------------------------------------------------\n")
            file_handle.write(f"# --- WINGLET TRANSITION START ---\n")
            
            wing_sweep_LE = np.deg2rad(wing_sweep_LE_deg)
            x_start, y_start, z_start = start_point
            x_lst = []
            y_lst = []
            z_lst = [] 

            #start to create the part that is "rotating", then only do one more for the "straigth part"
            #it also said that the c is cstart, so there cannot be any taper, and also the sweep is the one at the end...
            segment_length = np.linspace(0, transition_length, n_segments)  # so 1, 2, 3, 4, incremental
            
            dL = transition_length / n_segments
            #initialize everything as "normal" and then add the pieces
            x_station = x_start
            y_station = y_start
            z_station = z_start
            for i in range(1, n_segments) :
                #remember that all this is just a planform analysis, we don't actually move the wing
                #in x, only the sweep changes the position, since no taper
                # in y, only the delta transition changes the position
                # in z , what changes is the delta cant angle..
                cant_i_deg = i/n_segments * (90 - self.cant_angle)
                cant_i =  np.deg2rad(cant_i_deg)  #so 1/10, 2/10, etc
                #now move for the sweep 
                
            
                dy = dL * np.cos(cant_i)
                dx = dy * np.sin(wing_sweep_LE)  #the x is not on how much I move with dL but how much I actually move (dy)
                dz = dL * np.sin(cant_i)
                x_station -= dx  #every time I will take a piece away
            
                y_station += dy
                z_station += dz
                x_lst.append(x_station)
                y_lst.append(y_station)
                z_lst.append(z_station)
                
                #now we want to write automatically, but remember that the chord will NOT change
                file_handle.write('#\n')
                file_handle.write(f"#-------------------------------------------------------------\n")
                file_handle.write("SECTION\n")
                file_handle.write("#Xle    Yle    Zle     Chord   Ainc  Nspanwise  Sspace\n")
                file_handle.write(f"{x_station:.4f}  {y_station:.4f}  {z_station:.4f}   {self.root_c:.4f}   0.0   4   0\n")

            #now we want to create the last part of the winglet. it is only one more section (Straight)
            x_start_winglet = x_lst[-1]
            y_start_winglet = y_lst[-1]
            z_start_winglet = z_lst[-1]
            file_handle.write(f"#-------------------------------------------------------------\n")
            file_handle.write(f"# --- WINGLET STRAIGHT PART ---\n")
            file_handle.write(f"#-------------------------------------------------------------\n")

            
            x_end_winglet = x_start_winglet - self.b*np.sin(np.deg2rad(self.delta_le))
            y_end_winglet = y_start_winglet + self.b * np.sin(np.deg2rad(self.cant_angle))
            z_end_winglet = z_start_winglet + self.b * np.cos(np.deg2rad(self.cant_angle))

            #the ainc is specified the last line, which is the TWIST. 
            ainc = self.tip_twist - self.root_twist
            file_handle.write("SECTION\n")
            file_handle.write("#Xle    Yle    Zle     Chord   Ainc  Nspanwise  Sspace\n")
            file_handle.write(f"{x_end_winglet:.4f}  {y_end_winglet:.4f}  {z_end_winglet:.4f}   {self.tip_c:.4f}   {ainc}   20   0\n")

            print('writing on file completed')

                
        return 
    

    def build_geom_copy(self, n_segments, start_point, wing_sweep_LE_deg):  # I define a lot of segments to avoid the 
        #90 deg sharp turn (should not make a difference I believe (no interference turbulence, AVL) but slightly less Cl, less surface exp)
        transition_length =  0.2  #random val we will see
        winglet_file_path = r"C:\copiarefiles\aircraft_ady\ADY_AVL\winglet.txt"
    

        with open(winglet_file_path, 'w') as file_handle:
            file_handle.write(f"#-------------------------------------------------------------\n")
            file_handle.write(f"# --- WINGLET TRANSITION START ---\n")
            
            wing_sweep_LE = np.deg2rad(wing_sweep_LE_deg)
            x_start_winglet, y_start_winglet, z_start_winglet = start_point
            x_lst = []
            y_lst = []
            z_lst = [] 

            #start to create the part that is "rotating", then only do one more for the "straigth part"
            #it also said that the c is cstart, so there cannot be any taper, and also the sweep is the one at the end...
            segment_length = np.linspace(0, transition_length, n_segments)  # so 1, 2, 3, 4, incremental
            
            dL = transition_length / n_segments
            #initialize everything as "normal" and then add the pieces
            
            file_handle.write(f"#-------------------------------------------------------------\n")
            file_handle.write(f"# --- WINGLET STRAIGHT PART ---\n")
            file_handle.write(f"#-------------------------------------------------------------\n")

        
            x_end_winglet = x_start_winglet - self.b*np.sin(np.deg2rad(self.delta_le))
            y_end_winglet = y_start_winglet + self.b * np.sin(np.deg2rad(self.cant_angle))
            z_end_winglet = z_start_winglet + self.b * np.cos(np.deg2rad(self.cant_angle))

            #the ainc is specified the last line, which is the TWIST. 
            ainc = self.tip_twist - self.root_twist
            file_handle.write("SECTION\n")
            file_handle.write("#Xle    Yle    Zle     Chord   Ainc  Nspanwise  Sspace\n")
            file_handle.write(f"{x_end_winglet:.4f}  {y_end_winglet:.4f}  {z_end_winglet:.4f}   {self.tip_c:.4f}   {ainc}   20   0\n")

            print('writing on file completed')

                
        return 
    
    
    def build_wing_1w(self, n_segments, start_point, wing_sweep_LE_deg, wetted_area = False):
        #this function just appends the files

        old_reference_chars = '#Sref    Cref    Bref\n'
        old_reference_vals = [35,     1.809,     20.0]

        self.build_geom_copy(n_segments, start_point, wing_sweep_LE_deg)
        cant_angle = self.cant_angle

        wing_1w_path = fr"C:\copiarefiles\aircraft_ady\ADY_AVL\wing_1w_as4_wetted_area_{cant_angle}.txt"
        wing_1_path = r"C:\copiarefiles\aircraft_ady\ADY_AVL\wing_1_as4.txt"
        winglet_path = r"C:\copiarefiles\aircraft_ady\ADY_AVL\winglet.txt"

        files_to_merge = [wing_1_path, winglet_path
                          ]
        
        with open(wing_1w_path, 'w') as outfile:
            for file_path in files_to_merge:
                with open(file_path, 'r') as infile:
                    # Scrive il contenuto del file corrente nel file finale
                    outfile.write(infile.read())
                    
                    # Aggiunge una riga vuota tra i file per sicurezza (opzionale)
                    outfile.write("\n")

        if wetted_area:
            #read the file again and change the surface area to the projected area. we don't care about chord...
            added_b = 2 * self.b * np.sin(np.deg2rad(self.cant_angle))  # for example with a 0 can't angle there is 0 added b. ADD TWICE, specular !!!
            start_b = old_reference_vals[-1] #this is 20
            final_b = added_b + start_b

            added_area = (self.root_c + self.tip_c) * added_b /2       # (b1 + b2 )* h/2 , the chord is already double!
            start_area = old_reference_vals[0]
            final_area = added_area + start_area

            new_line_content = str(final_area) + '     ' + str(old_reference_vals[1]) + '     ' + str(final_b) + '\n'
            lines = []
            #now we update the file 
            with open(wing_1w_path, 'r') as f:
                lines = f.readlines()
            
            with open(wing_1w_path, 'w') as f:
                replace_next = False
                for i, line in enumerate(lines):
                    clean_line = line.strip().upper()

                    if replace_next:
                        f.write(new_line_content)
                        replace_next = False  #basically activate only once
                        continue

                    f.write(line)
                    # Se la riga inizia con la parola chiave (es. Sref)
                    if 'SREF' in clean_line and 'BREF' in clean_line:
                        # we need to skip this one , the one to modify is the next
                        replace_next = True
                




        _, wing_1w_name = os.path.split(wing_1w_path)
        return wing_1w_name


def run_AVL_analysis(avl_dir, input_geom, Cl, cant_angle, wetted_area = False):
    
    avl_exe = "avl352.exe"
    #define the commands to launch
    input_file_path = os.path.join(avl_dir, "temp_input.txt")
    if wetted_area: 
            input_commands = (
        f' load {input_geom}\n'
        f'oper\n'
        f'a c {Cl}\n'
        f'x\n'
        f'ft results_wing_1w_wetted_area_{cant_angle}.txt\n'
        f'quit\n'
        )
    else:
        input_commands = (
        f' load {input_geom}\n'
        f'oper\n'
        f'a c {Cl}\n'
        f'x\n'
        f'ft results_wing_1w_{cant_angle}.txt\n'
        f'quit\n'
        )

    # Execute AVL
    with open(input_file_path, "w") as f:
        f.write(input_commands)

    # Launch_process
    process = subprocess.run(f"{avl_exe} < temp_input.txt", shell=True, cwd = avl_dir, capture_output=True, text=True)

    print("--- AVL STDOUT ---")
    print(process.stdout)
    
    if process.stderr:
        print("--- AVL ERROR ---")
        print(process.stderr)

    if "File not found" in process.stdout or "Error" in process.stdout:
        print("ATTENZIONE: AVL ha riscontrato un errore nel caricamento o nel calcolo!")


def read_Cdi_vs_cant_angle(search_path):
    cd_ind_lst = []
    angle_lst = []
    for root, dir, filetupl in os.walk(search_path):
        for file in filetupl:
            if 'results' in file:
                full_path = os.path.join(root, file)
                with open(full_path, "r") as f:
                    text = f.read()

                #
                cdind_match = re.search(r"CDind\s*=\s*([0-9.+-Ee]+)", text)
                e_match     = re.search(r"e\s*=\s*([0-9.+-Ee]+)", text)

                cdind = float(cdind_match.group(1)) if cdind_match else None
                e     = float(e_match.group(1)) if e_match else None
                cd_ind_lst.append(cdind)
                angle_lst.append(file.split('_')[-1].split('.txt')[0])

    return angle_lst, cd_ind_lst


def read_Cdi_normalwing(search_path):
    cd_ind_lst = []
    angle_lst = []
    with open(search_path, "r") as f:
                    text = f.read()

                #
    cdind_match = re.search(r"CDind\s*=\s*([0-9.+-Ee]+)", text)
    e_match     = re.search(r"e\s*=\s*([0-9.+-Ee]+)", text)

    cdind = float(cdind_match.group(1)) if cdind_match else None
    e     = float(e_match.group(1)) if e_match else None

    return 0, cdind

                
    
if __name__ == '__main__':
    rag = np.deg2rad(12)
    print(np.tan(rag))
    wing_1 = Wing(b_1 = 5, b_2 = 10, root_c=2, tip_c=1, dihedral=0) 
    mac = wing_1.MAC()
    wing_LE_sweep = wing_1.calculate_sweep_deg()
    S_ref, x_ref, _, _ = wing_1.ady_centre_coords()
    avl_dir = r"C:\copiarefiles\aircraft_ady\ADY_AVL"
    CL_analysis = 0.6
    b_a320 = 14980   #this is transformed into 10 ms
    c_a320 = 1600
    b_winglet = 2106
    tip_c_winglet = 435  
    start_point_winglet = [0, 10, 0]
    #now I will normalize it   x_scaled : x_real = wing_1 : a320
    b_norm_winglet = b_winglet * 10 / b_a320
    c_root_norm_winglet = 1
    taper_winglet = tip_c_winglet / c_a320
    c_tip_norm_winglet = c_root_norm_winglet * taper_winglet
    le_angle = 67.1 
    cant_angle_lst = np.arange(0, 91, 5)
    twist_angle = -2

    path_wetted = r"C:\copiarefiles\aircraft_ady\ADY_AVL\results_wing1w_wetted_area"
    
    
    path_small = r"C:\copiarefiles\aircraft_ady\ADY_AVL\results_wing1w_smallaera"

    angles, cdinds = read_Cdi_vs_cant_angle(path_small)
    angles = [float(a) for a in angles]

    angles_sorted, cdinds_sorted = zip(*sorted(zip(angles, cdinds)))


    angle_0, cdi_0  = read_Cdi_normalwing(r"C:\copiarefiles\aircraft_ady\ADY_AVL\wing_normal_res")
    plt.figure(figsize=(8, 5), dpi=120)

    plt.plot(
        angles_sorted, cdinds_sorted,
        linestyle='-',
        color = 'green',
        linewidth=2.0,
        marker='s',
        markersize=6,
        label='Cdi_wing1w'
    )

    plt.plot(
        0, 0.0101876,
        linestyle='-',
        color = 'black',
        linewidth=2.0,
        marker='o',
        markersize=6,
        label='Cdi_wing2'
    )

    plt.plot(
        0, 0.0051380,
        linestyle='-',
        color = 'blue',
        linewidth=2.0,
        marker='o',
        markersize=6,
        label='Cdi_wing3'
    )
    plt.plot(
        0, cdi_0,
        linestyle='-',
        color = 'red',
        linewidth=2.0,
        marker='o',
        markersize=6,
        label='Cdi_wing1'
    )

    plt.plot()

    # -----------------------------
    # LABEL & TITOLO
    # -----------------------------
    plt.xlabel('cant angle [deg]', fontsize=12)
    plt.ylabel('CDi', fontsize=12)
    plt.title('Induced Cd with cant angle variation (same Sref) ', fontsize=14, fontweight='bold')

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

    # -----------------------------
    # SHOW / SAVE
    # -----------------------------
    plt.show()