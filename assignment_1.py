import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
from numpy import log as ln  #for reading simplicity
from mpl_toolkits.mplot3d import Axes3D
import sympy  #very useful to do derivatives and double check values
import matplotlib.path as mpath

#note, this script will use atmospheric pressure and 15 C for simplicity as standard

class PlotManager:
    def __init__(self, X, Y, x0, y0):
        """Inizializza con la griglia di calcolo comune"""
        self.X = X
        self.Y = Y
        self.x0 = x0
        self.y0 = y0

    def get_body_mask(self, psi_data, psi_0):
        """
        Ritorna una matrice booleana (Mask) dove True indica che il punto 
        è all'interno dell'ovale (psi = psi_0).
        """
        # 1. Estraiamo il path dal contour (senza plottare nulla)
        temp_fig, temp_ax = plt.subplots()
        cont = temp_ax.contour(self.X, self.Y, psi_data, levels=[psi_0])
        
        # Prendiamo il primo path trovato
        path = cont.get_paths()[0]
        plt.close(temp_fig)

        # 2. Creiamo una lista di tutti i punti della griglia (X, Y)
        # flattened_points sarà una matrice (1000000, 2)
        grid_points = np.column_stack((self.X.ravel(), self.Y.ravel()))

        # 3. Verifichiamo quali punti sono dentro il path
        mask_flat = path.contains_points(grid_points)
        
        # 4. Riportiamo la maschera alla forma della griglia (1000, 1000)
        return mask_flat.reshape(self.X.shape)

    def plot_phi_psi(self, phi, psi, title="Flow Field"):
        """Creates the plot for Phi e Psi"""
        data_to_plot = [phi, psi]
        labels = ['Potential $\\phi$', 'Stream Function $\\psi$']
        
        fig, axs = plt.subplots(1, 2, figsize=(12, 7))
        
        for i, ax in enumerate(axs):
            # creates contour with data from the function AS INPUT!!
            cs = ax.contour(self.X, self.Y, data_to_plot[i], 
                             levels=200, colors='cyan', linewidths=1)
            
            # personalization
            ax.set_title(f'{title} - {labels[i]}')
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.axis('equal')
            

        plt.tight_layout()
        plt.show()

    def pressure_and_velocity_with_arrows(self, plot_dict, R_cyl=None):
        u = plot_dict.get('Velocity x: u')
        v = plot_dict.get('Velocity y: v')
        vmag = plot_dict.get('Velocity Magnitude: vmag')
        cp = plot_dict.get('Pressure coefficient: cp')

        fig, axs = plt.subplots(1, 2, figsize=(12, 6), constrained_layout = True)
        axs = np.atleast_1d(axs)  # 🔒 rende axs sempre iterabile

        # ===== Disegno cilindro =====
        if R_cyl is not None:
            for ax in axs:
                circle = plt.Circle(
                    (self.x0, self.y0),
                    R_cyl,
                    color='white',
                    zorder=10
                )
                circle_b = plt.Circle(
                    (self.x0, self.y0),
                    R_cyl,
                    color='black',
                    fill=False,
                    lw=2,
                    zorder=11
                )
                ax.add_patch(circle)
                ax.add_patch(circle_b)

        # ===== VELOCITY FIELD =====
        levels = np.linspace(np.nanmin(vmag), np.nanmax(vmag), 15)
        cf0 = axs[0].contourf(self.X, self.Y, vmag, levels=levels, cmap='viridis')
        axs[0].contour(self.X, self.Y, vmag, levels=levels, colors='black', linewidths=0.15)
        axs[0].streamplot(
            self.X, self.Y,
            u / (vmag + 1e-12),
            v / (vmag + 1e-12),
            color='black',
            density=2, linewidth = 0.4
        )

        axs[0].set_title('normalized Velocity field')
        axs[0].set_xlabel('normalized space')
        axs[0].set_ylabel('normalized velocity magnitude [m/s]')
        axs[0].set_aspect('equal')
        fig.colorbar(cf0, ax = axs[0], shrink=0.6,
    pad=0.02,
    aspect=30)


        # ===== PRESSURE COEFFICIENT =====
        levels = np.linspace(np.nanmin(cp), np.nanmax(cp), 20)
        cf1 = axs[1].contourf(self.X, self.Y, cp, levels=levels, cmap='viridis', )
        axs[1].contour(self.X, self.Y, cp, levels=levels, colors='black', linewidths=0.5)

        axs[1].set_title('Pressure coefficient Cp')
        axs[1].set_xlabel('normalized space')
        axs[1].set_ylabel('Cp')
        axs[1].set_aspect('equal')
        fig.colorbar(cf1, ax = axs[1], shrink=0.6,
    pad=0.02,
    aspect=30)
        plt.show()


    def plot_grid(self, plots_dict,  detail=True, cols=3, figsize=(12, 10)):
        # 1. Recupero dati necessari
        psi_data = plots_dict.get('Psi')
        u_data = plots_dict.get('u')
        v_data = plots_dict.get('v')
        
        body_coords = None

        # 2. Identificazione della sagoma del corpo (una sola volta)
        if detail and psi_data is not None and u_data is not None:
            v_mag_sq = u_data**2 + v_data**2
            # Cerchiamo il ristagno
            idx = np.unravel_index(np.argmin(v_mag_sq), v_mag_sq.shape)
            psi_0 = psi_data[idx]

            # Estraiamo i punti della isoline senza creare maschere
            temp_fig, temp_ax = plt.subplots()
            cont = temp_ax.contour(self.X, self.Y, psi_data, levels=[psi_0])
            paths = cont.get_paths()
            if paths:
                # Prendiamo i punti (vertici) della linea di ristagno
                body_coords = paths[0].vertices 
            plt.close(temp_fig)

        # 3. Creazione della griglia di plot
        n_plots = len(plots_dict)
        rows = (n_plots + cols - 1) // cols
        fig, axs = plt.subplots(rows, cols, figsize=figsize)
        axes_flat = axs.flatten() if n_plots > 1 else [axs]

        for ax, (title, data) in zip(axes_flat, plots_dict.items()):
            # Plot del campo fluido (senza maschere complesse)
            cp = ax.contourf(self.X, self.Y, data, levels=10, cmap='viridis', alpha=0.8)
            ax.contour(self.X, self.Y, data, levels=50, colors='black', linewidths=0.5, alpha=0.3)
            
            # 4. IL "CEROTTO": Copriamo il corpo con un riempimento bianco
            if body_coords is not None:
                # Disegna il poligono bianco sopra i dati
                ax.fill(body_coords[:,0], body_coords[:,1], color='white', zorder=10)
                # Disegna il bordo del corpo più marcato
                ax.plot(body_coords[:,0], body_coords[:,1], color='black', linewidth=1.5, zorder=11)

            fig.colorbar(cp, ax=ax, shrink=0.8)
            ax.set_title(title)
            ax.set_aspect('equal')
            ax.set_xlim(np.min(self.X), np.max(self.X))
            ax.set_ylim(np.min(self.Y), np.max(self.Y))

        # Nascondi subplot vuoti
        for j in range(n_plots, len(axes_flat)):
            axes_flat[j].axis('off')

        plt.tight_layout()
        plt.show()


    def plot_multiple_plots(self, plots_dict, R_cyl=None):
        n_plots = len(plots_dict)
        cols = 3
        rows = (n_plots + cols - 1) // cols
        
        fig, axs = plt.subplots(rows, cols, figsize=(12, 6))
        axes_flat = axs.flatten() if n_plots > 1 else [axs]

        for ax, (title, data) in zip(axes_flat, plots_dict.items()):
            # Usiamo livelli automatici ma robusti
            # Se ci sono nan, contourf li gestisce bene (lascia bianco)
            cp = ax.contourf(self.X, self.Y, data, levels=25, cmap='viridis')
            fig.colorbar(cp, ax=ax)
            ax.contour(self.X, self.Y, data, levels=25, colors='black', linewidths=0.5, alpha=0.5)
            
            # Disegna il cerchio del cilindro
            if R_cyl is not None:
                circle = plt.Circle((self.x0, self.y0), R_cyl, color='white', zorder=10)
                ax.add_patch(circle)
                circle_b = plt.Circle((self.x0, self.y0), R_cyl, color='black', fill=False, lw=2, zorder=11)
                ax.add_patch(circle_b)

            ax.set_title(title)
            ax.set_aspect('equal')

        for j in range(n_plots, len(axes_flat)):
            axes_flat[j].axis('off')
        
        plt.tight_layout()
        plt.show()


class FlowElement:
    '''generic class that will keep count of general properties
    the first formulas are just placeholders. basically according to the theory of 
    OOP I should have them just in case. if I want to iterate on more instances, maybe
    I called the get_velocity in another way, and this will raise the error
    Instead, the get_local_pressure is the same for all, since the actual formula is the same for all
    this way I can plot the local pressure BOTH by creating a one element flowfield AND by
    using the local method.
    since the method needs the LOCAL velocities, I cannot simply say get_velocity because
    this would raise the method here, and basically I could not have different velocities calculations
    (in reality, in this case I would get NotimplementedError)'''


    def __init__(self):
        self.T0 = 15 + 273.15 # K 
        self.P0 = 101325 # Pa
        self.rho0 = 1.225 #standard air

    
    def get_velocity(self, X, Y):
        raise NotImplementedError
    
    def get_phi(self, X, Y):
        raise NotImplementedError

    def get_psi(self, X, Y):
        raise NotImplementedError
    
    def get_local_pressure_coeff(self, u, v, Uinf = 0):
        '''in this case I pass u, v and the undisturbed velocity Uinf. 
        in this case I can simply say patm + uinf ^2 = p + uloc^2 '''
        #  u, v = self.get_velocity DON'T DO!! if not it will call the method here,
        #raising a error. I need to pass the velocities as INPUTS!!
        V =np.sqrt( u**2 + v**2) # magnitude. v and u should be vectors. if at any points I do bernoulli should work

        #the local pressure is the atm + Uinf^2 - Uloc^2  (the faster the smaller correct)
        cp_loc = 1 - (V / Uinf)**2 if Uinf!=0 else 0
        return cp_loc  #since all is vectors, this should also be a matrix of scalars 


class Singularity(FlowElement):
    '''this class is the father for sources, sinks, vortices'''

    def __init__(self, strength, x0, y0):
        super().__init__()
        self.strength = strength
        self.x0 = x0
        self.y0 = y0
    
    def get_R(self, X, Y):
        '''this is good because I don't have to "calculate again" R in all classes,
        I can directly do "self.get_R" in both without rewriting it"
        '''
        x = X - self.x0
        y = Y - self.y0
        R = x**2 + y**2 + 1e-10 #smart, add a tiny value to avoid 0 (for plotting and ratios)
        return x, y, R   


class SourceSink(Singularity):
    def __init__(self, m, x0, y0):
        super().__init__(strength=m, x0=x0, y0=y0)
        self.m = m   #this is just kinda of an alias ,but makes it more readable for me
    def get_phi(self,X,Y):
        _, _, R = self.get_R(X,Y)
        return (self.m / (2 * np.pi) * ln(np.sqrt(R)))

    def get_psi(self, X, Y):
        #important, remember to use atan2 because works in all quadrants
        x, y, _ = self.get_R(X,Y)
        return  (self.m / (2 * np.pi) * np.arctan2(y, x))

    def get_velocity(self, X, Y):
        x,y,R = self.get_R(X,Y)
        u = (self.m / (2*np.pi)) / R * x
        v = (self.m / (2*np.pi)) / R * y
        return u, v 


class Vortex(Singularity):
    def __init__(self, Gamma, x0, y0):
        super().__init__(strength=Gamma, x0=x0, y0=y0)
        self.Gamma = Gamma

    def get_psi(self,X,Y):
        _,_, R = self.get_R(X,Y) 
        return (- self.Gamma / (2 * np.pi) * ln(np.sqrt(R)))

    def get_phi(self, X, Y):
        x,y,_ = self.get_R(X,Y)
        #important, remember to use atan2 because works in all quadrants
        return  (self.Gamma / (2 * np.pi) * np.arctan2(y, x))
    
    def get_velocity(self, X, Y):
        x,y,R = self.get_R(X,Y)
        u = (self.Gamma / (2*np.pi)) / R * y
        v = - (self.Gamma / (2*np.pi)) / R * x
        return u, v 


class FreeFlow(FlowElement):
    def __init__(self, alpha, U):
        super().__init__()
        self.alpha = alpha
        self.U = U
    
    def deg_to_rad(self):
        return  np.radians(self.alpha) 
    
    def get_phi(self, X, Y):
        alpha_rad = self.deg_to_rad()
        return  (self.U * (X*np.cos(alpha_rad) + Y* np.sin(alpha_rad)))
    
    def get_psi(self, X, Y):
        alpha_rad = self.deg_to_rad()
        return (self.U * (Y*np.cos(alpha_rad) - X* np.sin(alpha_rad)))

    def get_velocity(self, X, Y):
        alpha_rad = self.deg_to_rad()

        #careful not to return a scalar value
        u = self.U * np.cos(alpha_rad) * np.ones_like(X)
        v = self.U * np.sin(alpha_rad)  * np.ones_like(X)
        return u, v


class FlowSystem:

    '''this is the class that will do all the calculations on the code, adding the phi
    and the psi, finding the overall velocity field'''  
    def __init__(self, elemental_flows_lst, x, y):
        self.elemental_flows_lst = elemental_flows_lst
        self.X, self.Y = np.meshgrid(x,y)

    def get_velocity(self,X,Y):
        ''' it is important to do elementwise sum, bugt python already does this normally
        this means that velocity + velocity 2 will act on the elements'''
        u_tot = np.zeros_like(X)  #this is important to avoid concat. it is because we use numpy
        v_tot = np.zeros_like(X) #do not have to use y, they are both same
        for el in self.elemental_flows_lst:
            u, v = el.get_velocity(X,Y)
            u_tot += u
            v_tot += v
        return u_tot, v_tot

    def get_pressure_coeff(self, X,Y, Uinf=0):
        P0 = self.elemental_flows_lst[0].P0  #get it just once, pressure is the SAME
        rho0 = self.elemental_flows_lst[0].rho0
        u_tot, v_tot = self.get_velocity(X,Y) 
        V = np.sqrt(u_tot**2 + v_tot**2)

        #using a system of eqn for bernoulli we can avoid using the rho
        cp = 1 - (V / Uinf)**2 if Uinf!= 0 else 0
        return cp
    
    def get_phi(self,X,Y):
        phi_tot = np.zeros_like(X)
        for el in self.elemental_flows_lst:
            phi = el.get_phi(X,Y)
            phi_tot += phi
        return phi_tot
    
    def get_psi(self,X,Y):
        psi_tot = np.zeros_like(X)
        for el in self.elemental_flows_lst:
            psi = el.get_psi(X,Y)
            psi_tot += psi
        return psi_tot
    def get_all_functions(self, X,Y, Uinf = 0):
        phi = self.get_phi(X, Y)
        psi = self.get_psi(X,Y)
        u, v = self.get_velocity(X,Y)
        press = self.get_pressure_coeff(X,Y, Uinf)
        vmag = np.sqrt(u**2 + v**2)
        return phi, psi, u, v, vmag, press

    
    def get_value_at(self, x_user, y_user, data_matrix, x_vec, y_vec):
        """ 
        Finds the closest point in the matrix to what the user passed.
        """
        # Trova l'indice più vicino per x e y
        idx_x = np.abs(x_vec - x_user).argmin()
        idx_y = np.abs(y_vec - y_user).argmin()
        
        # Restituisce il valore (nota: nelle matrici l'asse y è la riga, x la colonna)
        return data_matrix[idx_y, idx_x]
    
    def add_element(self, element):
        self.elemental_flows_lst.append(element)


class Doublet(Singularity):
    def __init__(self, kappa, x0, y0):
        super().__init__(strength=kappa, x0=x0, y0=y0)
        self.kappa = kappa

    def get_psi(self, X, Y):
        x, y, R = self.get_R(X, Y)
        # Psi = - (kappa / 2pi) * (y / (x^2 + y^2))
        return - (self.kappa / (2 * np.pi)) * (y / R)

    def get_phi(self, X, Y):
        x, y, R = self.get_R(X, Y)
        # Phi = (kappa / 2pi) * (x / (x^2 + y^2))
        return (self.kappa / (2 * np.pi)) * (x / R)

    def get_velocity(self, X, Y):
        x, y, R = self.get_R(X, Y)
        # Derivate analitiche del dipolo
        factor = self.kappa / (2 * np.pi)
        u = - factor * (x**2 - y**2) / (R**2 + 1e-10) # Corretto segno per dipolo controcorrente
        v = - factor * (2 * x * y) / (R**2 + 1e-10)
        return u, v    


class RankineOval(FlowSystem):
    def __init__(self, m, x0_dist, y0, Uinf, alpha, x, y):
        source = SourceSink(m = m, x0 = -x0_dist, y0=y0)
        sink = SourceSink(m = m, x0 = x0_dist, y0=y0)
        freestream = FreeFlow(alpha = alpha, U=Uinf)

        #get the properties of system
        super().__init__([source, sink, freestream], x, y)


class CylinderFlow(FlowSystem):
    def __init__(self, Uinf, alpha, R_cyl, Gamma, x, y, x0, y0):
        # 1. Teorema del Cerchio: Kappa = 2 * pi * Uinf * R^2
        # Questo assicura che psi=0 sulla superficie del cilindro
        kappa = 2 * np.pi * Uinf * (R_cyl**2)
        self.x0 = x0
        self.y0 = y0
        doublet = Doublet(kappa=kappa, x0=x0, y0=y0)
        freestream = FreeFlow(alpha=alpha, U=Uinf)
        
        elements = [freestream, doublet]
        
        if Gamma != 0:
            # Vortice orario positivo o negativo a seconda della convenzione
            # Qui assumiamo Gamma positivo = circolazione standard
            vortex = Vortex(Gamma=Gamma, x0=x0, y0=y0)
            elements.append(vortex)
            
        super().__init__(elements, x, y)
        self.R_cyl = R_cyl # Salviamo il raggio per il masking

    def mask_inside_cylinder(self, X, Y, data_lst, shrink=0.98, ):

        mask =  (X - self.x0)**2 + (Y - self.y0)**2 < (self.R_cyl * shrink)**2
        for data in data_lst:
            data[mask] = np.nan
        return data_lst
    def get_lift_coeff(self):
        rho0 = self.elemental_flows_lst[0].rho0
        Lift = rho0*Gamma*Uinf
        return  Lift / (0.5 * rho0 * self.R_cyl * 2 * Uinf**2)


def adimensionalisation(R_adi, Lx=600.0, Ly=400.0, Uinf=None, x0=None, y0=None, m=None, Gamma=None):
    # dimensionless reference

    L_ref = max(Lx, Ly)
    if x0 > Lx or y0>Ly:
        print('you are placing elements outside the bounds!, the domain is 300x200, put smaller coordinates!')
        return
    Uinf_adi = 1.0 if Uinf is not None else None

    x0_adi = x0 / L_ref if x0 is not None else None
    y0_adi = y0 / L_ref if y0 is not None else None

    m_adi = m / (Uinf * L_ref) if all(v is not None for v in [m, Uinf, L_ref]) else None
    Gamma_adi = Gamma / (Uinf * 2*R_adi) if all(v is not None for v in [Gamma, Uinf, L_ref]) else None

    # grid in dimensionless box
    x = np.linspace(-Lx/2, Lx/2, 500) / L_ref
    y = np.linspace(-Ly/2, Ly/2, 500) / L_ref

    X, Y = np.meshgrid(x, y)

    return x0_adi, y0_adi, x, y, X, Y, m_adi, Uinf_adi, Gamma_adi


def leveled_flight(m, Uinf, R_cyl, b_semi):
    '''this is for the magnus effect cylinder, automatically put at x = 0, y = 0'''
    
    rho0 = 1.225
    b = 2* b_semi
   
    m_cyl_tot = 2 * 100 * (np.pi * R_cyl**2 * b_semi)   #assume smaller cylinders, 100 kg per meter length and diameter
    d_cyl = 2*R_cyl
    

    W = (m ) * 9.81

    Uinf = 20
    x = np.linspace(-10, 10, 500 ) #plot downstream "wake" 10 times cyl
    y = np.linspace(-10, 10, 500)
    X,Y = np.meshgrid(x,y)
    #now we go backwards to define Gamma = L / (rho*Uinf),  L  * L_correction= w 
    k = 0.2   #since we know that with a aspect ratio of 5 it is "almost perfect", we consider almost perfect with 5% deviation. consider changing to other formula
    AR = b**2 / (b * d_cyl )   #very approximate estimate of the aspect ratio
    L_correction = 1 - k/AR

    L_per_span = W / b
    Gamma_tot = (L_per_span) / (rho0 * Uinf * L_correction)   #this is the circulation TOTAL (not L')

    Gamma_cyl = Gamma_tot / 2
    Omega = Gamma_cyl  / (2*np.pi * R_cyl**2) #find the rotating speed imposing the circulation

    print(f'rotational speed required to sustain flight = {Omega} rad/s')
    Lambda = Omega* R_cyl / Uinf

    C_torque = 0.02 if Lambda <2 else 0.05 if Lambda <5 else 0.09
    mech_eta = 0.9
    A_proj = np.pi * R_cyl**2
    Torque_aero = 0.5 * rho0 * Uinf**2 * A_proj * R_cyl * C_torque  # [N m]
    Power_aero = Torque_aero * Omega   # [W] - potenza che la rotazione deve bilanciare
    # Potenza elettrica (input) considerando efficienza meccanica:
    Power_elec = Power_aero / mech_eta
    print(f'Power required to sustain flight = {Power_aero} W')
    aircraft = CylinderFlow(Uinf, 0, R_cyl=R_cyl, Gamma=1, x = x, y = y, x0 = x0_adi, y0 = y0_adi )

    phi, psi, u, v, vmag, cp = aircraft.get_all_functions(X,Y,Uinf_adi)

    data_lst = [phi, psi, u, v, vmag, cp]
    [phi, psi, u, v, vmag, cp] = aircraft.mask_inside_cylinder(X, Y, data_lst)

    pm = PlotManager(X,Y, x0_adi, y0_adi)
    plots_dict = {
        'Potential: Phi': phi,
        'Stream Function: psi': psi,
        'Velocity Magnitude: vmag': vmag,
        'Pressure coefficient: cp': cp,
        'Velocity x: u' :u,
        'Velocity y: v' : v
    }
   # pm.plot_multiple_plots(plots_dict)
    return b, d_cyl, Omega, Power_aero

if __name__ == '__main__':
    #all quantities are now in a 300 x 200 box. DO NOT EXCEED +-300 for x and +-400 for y


    


#%%
    x0 = 0 
    y0 = 0
    m = 10
    Uinf = 10
    Gamma = 0.3
    alpha = 0
    R_adi = 0.05
    x0_adi, y0_adi, x, y, X, Y, m_adi, Uinf_adi, Gamma_adi = adimensionalisation(x0 =x0,y0 = y0, m= m,Uinf= Uinf,Gamma = Gamma, R_adi = R_adi)
    
    #all quantities are now in a 300 x 200 box. DO NOT EXCEED 1 with radius

    cyl_sys = CylinderFlow(Uinf=Uinf_adi, alpha = alpha, R_cyl= R_adi, Gamma= Gamma_adi, x=x, y=y, x0 = x0_adi, y0 = y0_adi)

    phi, psi, u, v, vmag, cp = cyl_sys.get_all_functions(X,Y,Uinf_adi)

    data_lst = [phi, psi, u, v, vmag, cp]
    plot_dict = {
        'Potential: Phi': phi,
        'Stream Function: psi': psi,
        'Velocity Magnitude: vmag': vmag,
        'Pressure coefficient: cp': cp,
        'Velocity x: u' :u,
        'Velocity y: v' : v
    }
    [phi, psi, u, v, vmag, cp] = cyl_sys.mask_inside_cylinder(X, Y, data_lst)


    Cl = cyl_sys.get_lift_coeff()
    print(Cl)

    # 6. Plotting
    pm = PlotManager(X, Y, x0_adi, y0_adi)
    
    plot_dict = {
        'Potential: Phi': phi,
        'Stream Function: psi': psi,
        'Velocity Magnitude: vmag': vmag,
        'Pressure coefficient: cp': cp,
        'Velocity x: u' :u,
        'Velocity y: v' : v
    }

    pm.plot_multiple_plots(plot_dict, R_adi)
    # Passiamo R_adi al plotter per disegnare il cerchio bianco pulito
    #pm.plot_multiple_plots(plot_dict, R_cyl=R_adi)

    #pm.pressure_and_velocity_with_arrows( plot_dict, R_cyl=R_adi)
#%%
    semispans = np.arange(0.2, 1, 0.1)
    omega_lst = []
    power_lst = []
    for semispan in semispans:
        b, r, omega, power = leveled_flight(m = 0.5, Uinf=20, R_cyl= semispan, b_semi = 0.5)
        omega_lst.append(power)
        power_lst.append(r)
    aa = []
    bb = []
    for semispan in semispans:
        b, r, omega, power = leveled_flight(m = 0.5, Uinf=20, R_cyl= 0.5, b_semi = semispan)
        print(power)
        aa.append(power)
        bb.append(b)


    plt.figure(figsize=(8, 5), dpi=120)

    plt.plot(
        bb, aa,
        linestyle='-',
        color = 'green',
        linewidth=2.0,
        marker='s',
        markersize=6,
        label='Power scaling span'
    )
    plt.plot(
        power_lst, omega_lst,
        linestyle='-',
        color = 'red',
        linewidth=2.0,
        marker='s',
        markersize=6,
        label='Power scaling diameter'
    )

    plt.plot()

    # -----------------------------
    # LABEL & TITOLO
    # -----------------------------
    plt.xlabel('span/diameter', fontsize=12)
    plt.ylabel('Power required', fontsize=12)
    plt.title('Variation of power required with span / diameter ', fontsize=14, fontweight='bold')

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