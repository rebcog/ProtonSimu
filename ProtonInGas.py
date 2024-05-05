#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import csv

colors = [ '#000080', '#0000cc', '#0033ff', '#1e90ff', '#4169e1']

# Constants GOOD
me = 0.511                              # MeV/c^2, electron mass
c = 3e10                                # cm/s, speed of light
K = 0.307075                            # MeV cm^2 / mol, constant in Bethe-Bloch formula
m_proton = 938.27                       # MeV/c^2, proton mass
z = 1                                    # charge of proton

#Water Constants GOOD
Z_water = 7.42                           # Effective atomic number for water (The physics of proton therapy Wayne D Newhauser and Rui Zhang))
A_water = 18.015                        # g/mol, molar mass of water
I_water = 75e-6                         # MeV, mean excitation energy for water
rho_water = 1.0                         # g/cm³, density of water

#Beam Constants
beam_radius = 0.5                                # radius of the proton beam in cm
cross_sectional_area = np.pi * beam_radius**2     # area in cm^2

#Math methods 

# Bethe-Bloch energy loss function
def bethe_bloch_energy(E_kin):
    beta = np.sqrt(1 - (m_proton / (E_kin + m_proton))**2)
    gamma = 1 / np.sqrt(1 - beta**2)
    Tmax = (2 * me * c**2 * beta**2 * gamma**2) /            (1 + (2 * gamma * me / m_proton) + (me / m_proton)**2)
    dEdx = K * Z_water * (z**2 / (beta**2 * A_water)) *            (np.log(2 * me * c**2 * beta**2 * gamma**2 * Tmax / I_water) - beta**2)
    return dEdx*rho_water                                         

def simulate_proton_beam(E, num_protons=2000, target_size=20, beam_offset=10): 
    depths_energy_deposition = []
    E_threshold = 0.001  # Energy threshold in MeV below which the simulation stops

    
    for _ in range(num_protons):
        depth = beam_offset
        E_current = E
        while E_current > E_threshold and depth < (beam_offset + target_size):
            dEdx = bethe_bloch_energy(E_current)
            mean_step_length = 1 / dEdx
            step_length = -np.log(np.random.uniform()) * mean_step_length
            mean_delta_E = dEdx * step_length
            delta_E = np.random.normal(mean_delta_E, .05 * mean_delta_E)

            depth += step_length
            if beam_offset <= depth <= (beam_offset + target_size):
                # Convert delta_E from MeV/cm to MeV/cm^3 by dividing by the cross-sectional area
                energy_density = delta_E / cross_sectional_area
                depths_energy_deposition.append((depth - beam_offset, energy_density))

            E_current -= delta_E
            if E_current < 0:
                E_current = 0

    return depths_energy_deposition

# Initialize plot
plt.figure(figsize=(10, 6))


# Simulate and plot for each energy level
energies = np.arange(50, 251, 50)  
for i, E in enumerate(energies):
    result = simulate_proton_beam(E)
    depths, energy_deposits = zip(*result)
    depth_bins = np.linspace(0, max(depths), 100)
    dose_profile, _ = np.histogram(depths, bins=depth_bins, weights=energy_deposits)
    plt.plot(depth_bins[:-1], dose_profile, label=f'{E} MeV', color=colors[i])  # Use color from list


#edit the graph to look bette
plt.xlabel('Penetration Depth (cm)', fontsize=13, labelpad=10)
plt.ylabel('Energy Deposition (MeV/cm³)', fontsize=13, labelpad=10)
plt.title('Proton Beam Energy Deposition in Water', fontsize=14, pad=20)
plt.legend(title="Proton Energy", fontsize='large', title_fontsize='large')
plt.tick_params(axis='both', which='major', labelsize=12)

# Use scientific notation for y-axis
ax = plt.gca()                                                 # Get current axis
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
y_formatter = ScalarFormatter(useMathText=True) 
y_formatter.set_scientific(True) 
y_formatter.set_powerlimits((-3,3))
ax.yaxis.set_major_formatter(y_formatter)

plt.grid(True)
plt.tight_layout(pad=2.0)                                      # Adjust the padding between and around subplots
plt.show()
#Saving file
# Constants and functions (insert your existing constants and functions here)

def save_data_to_csv(data, filename):
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Depth (cm)', 'Energy Deposition (MeV/cm^3)'])  # Writing headers
        writer.writerows(data)

# Simulate and save data for each energy level
energies = np.arange(50, 251, 50)  # Adjusted energy levels if necessary
for i, E in enumerate(energies):
    result = simulate_proton_beam(E)
    if result:  # Check if result is not empty
        depths, energy_deposits = zip(*result)
        depth_bins = np.linspace(0, max(depths), 100)
        dose_profile, _ = np.histogram(depths, bins=depth_bins, weights=energy_deposits)
        
        # Preparing data to save. Pair each depth bin with the corresponding dose profile value
        data_to_save = list(zip(depth_bins[:-1], dose_profile))
        
        # Save data to CSV, each energy will have a different file
        save_data_to_csv(data_to_save, f'/Users/reb/Desktop/Senior Thesis/proton_beam_{E}_MeV_energy_deposition.csv')
        


# In[17]:


#Bone Profile 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import csv

water_colors = [ '#000080', '#0000cc', '#0033ff', '#1e90ff', '#4169e1']
colors = ['#a53606', '#b32db5', '#881a58', '#0e288e', '#164c64']


# Constants
me = 0.511                           # MeV/c^2, electron mass
c = 3e10                             # cm/s, speed of light
K = 0.307075                         # MeV cm^2 / mol, constant in Bethe-Bloch formula
m_proton = 938.27                    # MeV/c^2, proton mass
z = 1                                # charge of proton

#Bone  #Good 
Z_bone = 13.2517                     # Effective atomic number for bone (cortical) 
A_bone = 21.332332406                # g/mol, molar mass of bone (https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=120)
I_bone = 1.064e-4                    # MeV, mean excitation energy for bone (https://physics.nist.gov/cgi-bin/Star/compos.pl?refer=ap&matno=120)
rho_bone = 1.85                      # g/cm³, density of Bone (https://physics.nist.gov/cgi-bin/Star/compos.pl?refer=ap&matno=120)

#Beam Constants
beam_radius = 0.5                                # radius of the proton beam in cm
cross_sectional_area = np.pi * beam_radius**2     # area in cm^2

#Math methods 

# Bethe-Bloch energy loss function
def bethe_bloch_energy(E_kin):
    beta = np.sqrt(1 - (m_proton / (E_kin + m_proton))**2)
    gamma = 1 / np.sqrt(1 - beta**2)
    Tmax = (2 * me * c**2 * beta**2 * gamma**2) /            (1 + (2 * gamma * me / m_proton) + (me / m_proton)**2)
    dEdx = K * Z_bone * (z**2 / (beta**2 * A_bone)) *            (np.log(2 * me * c**2 * beta**2 * gamma**2 * Tmax / I_bone) - beta**2)
    return dEdx*rho_bone                                         

def simulate_proton_beam(E, num_protons=2000, target_size=20, beam_offset=10): 
    depths_energy_deposition = []
    E_threshold = 0.001  # Energy threshold in MeV below which the simulation stops

    
    for _ in range(num_protons):
        depth = beam_offset
        E_current = E
        while E_current > E_threshold and depth < (beam_offset + target_size):
            dEdx = bethe_bloch_energy(E_current)
            mean_step_length = 1 / dEdx
            step_length = -np.log(np.random.uniform()) * mean_step_length
            mean_delta_E = dEdx * step_length
            delta_E = np.random.normal(mean_delta_E, .05 * mean_delta_E)

            depth += step_length
            if beam_offset <= depth <= (beam_offset + target_size):
                # Convert delta_E from MeV/cm to MeV/cm^3 by dividing by the cross-sectional area
                energy_density = delta_E / cross_sectional_area
                depths_energy_deposition.append((depth - beam_offset, energy_density))

            E_current -= delta_E
            if E_current < 0:
                E_current = 0

    return depths_energy_deposition

#NEW
def read_csv_data(filename):
    depths = []
    energy_depositions = []
    with open(filename, 'r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip the header row
        for row in csv_reader:
            depths.append(float(row[0]))
            energy_depositions.append(float(row[1]))
    return depths, energy_depositions
#END

# Initialize plot
plt.figure(figsize=(12, 7))

#NEW
energies = np.arange(50, 251, 50)
for i, E in enumerate(energies):
    filename = f'proton_beam_{E}_MeV_energy_deposition.csv'
    depths, energy_deposits = read_csv_data(filename)
    plt.plot(depths, energy_deposits, label=f'W {E} MeV', color=water_colors[i], alpha=0.3)
#END


#Bone Graph
# Simulate and plot for each energy level
energies = np.arange(50, 251, 50)  
for i, E in enumerate(energies):
    result = simulate_proton_beam(E)
    depths, energy_deposits = zip(*result)
    depth_bins = np.linspace(0, max(depths), 100)
    dose_profile, _ = np.histogram(depths, bins=depth_bins, weights=energy_deposits)
    plt.plot(depth_bins[:-1], dose_profile, label=f'B {E} MeV', color=colors[i])  # Use color from list
    



#edit the graph to look bette
plt.xlabel('Penetration Depth (cm)', fontsize=13, labelpad=10)
plt.ylabel('Energy Deposition (MeV/cm³)', fontsize=13, labelpad=10)
plt.title('Comparative Proton Beam Energy Deposition in Cortical Bone and Water', fontsize=14, pad=20)
plt.legend(title="Proton Energy", title_fontsize='large', ncol=2, fontsize='medium')
plt.tick_params(axis='both', which='major', labelsize=12)

# Use scientific notation for y-axis
ax = plt.gca()  # Get current axis
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
y_formatter = ScalarFormatter(useMathText=True) 
y_formatter.set_scientific(True) 
y_formatter.set_powerlimits((-3,3))
ax.yaxis.set_major_formatter(y_formatter)

plt.grid(True)
plt.tight_layout(pad=2.0)  # Adjust the padding between and around subplots
plt.show()


# In[18]:


#Lung Profile 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import csv

water_colors = [ '#000080', '#0000cc', '#0033ff', '#1e90ff', '#4169e1']

 
# Constants
me = 0.511                           # MeV/c^2, electron mass
c = 3e10                             # cm/s, speed of light
K = 0.307075                         # MeV cm^2 / mol, constant in Bethe-Bloch formula
m_proton = 938.27                    # MeV/c^2, proton mass
z = 1                                # charge of proton

#Lung Constants GOOD (https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=190)

Z_lung = 7.770249877812358                 # Effective atomic number for lung
A_lung = 14.196677382                      # g/mol, molar mass of lung (same)
I_lung = 75.3e-6                           # MeV, mean excitation energy for lung 
rho_lung = 1.05                            # g/cm³, density of lung 


#Beam Constants
beam_radius = 0.5                                # radius of the proton beam in cm
cross_sectional_area = np.pi * beam_radius**2     # area in cm^2

#Math methods 

# Bethe-Bloch energy loss function
def bethe_bloch_energy(E_kin):
    beta = np.sqrt(1 - (m_proton / (E_kin + m_proton))**2)
    gamma = 1 / np.sqrt(1 - beta**2)
    Tmax = (2 * me * c**2 * beta**2 * gamma**2) /            (1 + (2 * gamma * me / m_proton) + (me / m_proton)**2)
    dEdx = K * Z_lung * (z**2 / (beta**2 * A_lung)) *            (np.log(2 * me * c**2 * beta**2 * gamma**2 * Tmax / I_lung) - beta**2)
    return dEdx*rho_lung                                         

def simulate_proton_beam(E, num_protons=2000, target_size=20, beam_offset=10): 
    depths_energy_deposition = []
    E_threshold = 0.001  # Energy threshold in MeV below which the simulation stops

    
    for _ in range(num_protons):
        depth = beam_offset
        E_current = E
        while E_current > E_threshold and depth < (beam_offset + target_size):
            dEdx = bethe_bloch_energy(E_current)
            mean_step_length = 1 / dEdx
            step_length = -np.log(np.random.uniform()) * mean_step_length
            mean_delta_E = dEdx * step_length
            delta_E = np.random.normal(mean_delta_E, .05 * mean_delta_E)

            depth += step_length
            if beam_offset <= depth <= (beam_offset + target_size):
                # Convert delta_E from MeV/cm to MeV/cm^3 by dividing by the cross-sectional area
                energy_density = delta_E / cross_sectional_area
                depths_energy_deposition.append((depth - beam_offset, energy_density))

            E_current -= delta_E
            if E_current < 0:
                E_current = 0

    return depths_energy_deposition

#NEW
def read_csv_data(filename):
    depths = []
    energy_depositions = []
    with open(filename, 'r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip the header row
        for row in csv_reader:
            depths.append(float(row[0]))
            energy_depositions.append(float(row[1]))
    return depths, energy_depositions
#END

# Initialize plot
plt.figure(figsize=(12, 7))

#NEW
energies = np.arange(50, 251, 50)
for i, E in enumerate(energies):
    filename = f'proton_beam_{E}_MeV_energy_deposition.csv'
    depths, energy_deposits = read_csv_data(filename)
    plt.plot(depths, energy_deposits, label=f'W {E} MeV', color=water_colors[i], alpha=0.3)
#END

# Simulate and plot for each energy level
energies = np.arange(50, 251, 50)  
for E in energies:
    result = simulate_proton_beam(E)
    depths, energy_deposits = zip(*result)
    depth_bins = np.linspace(0, max(depths), 100)
    dose_profile, _ = np.histogram(depths, bins=depth_bins, weights=energy_deposits)
    plt.plot(depth_bins[:-1], dose_profile, label=f'L {E} MeV')

#edit the graph to look bette
plt.xlabel('Penetration Depth (cm)', fontsize=13, labelpad=10)
plt.ylabel('Energy Deposition (MeV/cm³)', fontsize=13, labelpad=10)
plt.title('Comparative Proton Beam Energy Deposition in Lung Tissue and Water', fontsize=14, pad=20)
plt.legend(title="Proton Energy", title_fontsize='large', ncol=2, fontsize='medium')
plt.tick_params(axis='both', which='major', labelsize=12)

# Use scientific notation for y-axis
ax = plt.gca()  # Get current axis
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
y_formatter = ScalarFormatter(useMathText=True) 
y_formatter.set_scientific(True) 
y_formatter.set_powerlimits((-3,3))
ax.yaxis.set_major_formatter(y_formatter)

plt.grid(True)
plt.tight_layout(pad=2.0)  # Adjust the padding between and around subplots
plt.show()


# In[62]:


## Calculate the molar mass of bone A_bone

atomic_masses = {
    'H': 1.008,
    'C': 12.011,
    'N': 14.007,
    'O': 15.999,
    'Mg': 24.305,
    'P': 30.974,
    'S': 32.065,
    'Ca': 40.078,
    'Zn': 65.38
}

weight_fractions = {
    'H': 0.047234,
    'C': 0.144330,
    'N': 0.041990,
    'O': 0.446096,
    'Mg': 0.002200,
    'P': 0.104970,
    'S': 0.003150,
    'Ca': 0.209930,
    'Zn': 0.000100
}

molar_mass_bone = sum(weight_fractions[element] * atomic_masses[element] for element in weight_fractions)
molar_mass_bone


# In[63]:


#Calculating effective atomic number of cortical (dense, low porosity)(Z)

atomic_numbers = {
    'H': 1,    # Hydrogen
    'C': 6,    # Carbon
    'N': 7,    # Nitrogen
    'O': 8,    # Oxygen
    'Mg': 12,  # Magnesium
    'P': 15,   # Phosphorus
    'S': 16,   # Sulfur
    'Ca': 20,  # Calcium
    'Zn': 30   # Zinc
}

weight_fractions = {
    'H': 0.047234,
    'C': 0.144330,
    'N': 0.041990,
    'O': 0.446096,
    'Mg': 0.002200,
    'P': 0.104970,
    'S': 0.003150,
    'Ca': 0.209930,
    'Zn': 0.000100
}

#commonly used is p = 3 for photoelectric interactions
p = 3

# Calculate the effective atomic number using the weight fractions and atomic numbers
Z_eff = (sum(weight_fractions[element] * (atomic_numbers[element] ** p) for element in weight_fractions)) ** (1/p)

print(Z_eff)


# In[66]:


#Molar mass of lung (A)


atomic_masses = {
    1: 1.008,    # Hydrogen
    6: 12.011,   # Carbon
    7: 14.007,   # Nitrogen
    8: 15.999,   # Oxygen
    11: 22.990,  # Sodium
    12: 24.305,  # Magnesium
    15: 30.974,  # Phosphorus
    16: 32.065,  # Sulfur
    17: 35.453,  # Chlorine
    19: 39.098,  # Potassium
    20: 40.078,  # Calcium
    26: 55.845,  # Iron
    30: 65.38    # Zinc
}

weight_fractions = {
    1: 0.101278,
    6: 0.102310,
    7: 0.028650,
    8: 0.757072,
    11: 0.001840,
    12: 0.000730,
    15: 0.000800,
    16: 0.002250,
    17: 0.002660,
    19: 0.001940,
    20: 0.000090,
    26: 0.000370,
    30: 0.000010
}

molar_mass_lung = sum(weight_fractions[atom] * atomic_masses[atom] for atom in weight_fractions)
molar_mass_lung


# In[67]:


#Effective Atomic Number for Lung (Z)

atomic_numbers = {
    1: 1,     # Hydrogen
    6: 6,     # Carbon
    7: 7,     # Nitrogen
    8: 8,     # Oxygen
    11: 11,   # Sodium
    12: 12,   # Magnesium
    15: 15,   # Phosphorus
    16: 16,   # Sulfur
    17: 17,   # Chlorine
    19: 19,   # Potassium
    20: 20,   # Calcium
    26: 26,   # Iron
    30: 30    # Zinc
}

weight_fractions = {
    1: 0.101278,
    6: 0.102310,
    7: 0.028650,
    8: 0.757072,
    11: 0.001840,
    12: 0.000730,
    15: 0.000800,
    16: 0.002250,
    17: 0.002660,
    19: 0.001940,
    20: 0.000090,
    26: 0.000370,
    30: 0.000010
}

#common choice is p=3 for photoelectric effect
p = 3

Z_eff = (sum(weight_fractions[Z] * (atomic_numbers[Z] ** p) for Z in weight_fractions)) ** (1/p)

Z_eff


# In[19]:


#Brain Profile

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import csv

water_colors = [ '#000080', '#0000cc', '#0033ff', '#1e90ff', '#4169e1']

colors = ['#980000', '#984500', '#a77f03', '#0c8900', '#004b83']
 
# Constants
me = 0.511                           # MeV/c^2, electron mass
c = 3e10                             # cm/s, speed of light
K = 0.307075                         # MeV cm^2 / mol, constant in Bethe-Bloch formula
m_proton = 938.27                    # MeV/c^2, proton mass
z = 1                                # charge of proton

#Brain Constants #https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=123

Z_brain =  7.75                # Effective atomic number for brain
A_brain = 14.03             # g/mol, molar mass of brain
I_brain =  73.3e-6           # MeV, mean excitation energy for brain
rho_brain = 1.03             # g/cm³, density of brain 


#Beam Constants
beam_radius = 0.5                                # radius of the proton beam in cm (5 mm) 
cross_sectional_area = np.pi * beam_radius**2     # area in cm^2

#Math methods 

# Bethe-Bloch energy loss function
def bethe_bloch_energy(E_kin):
    beta = np.sqrt(1 - (m_proton / (E_kin + m_proton))**2)
    gamma = 1 / np.sqrt(1 - beta**2)
    Tmax = (2 * me * c**2 * beta**2 * gamma**2) /            (1 + (2 * gamma * me / m_proton) + (me / m_proton)**2)
    dEdx = K * Z_brain * (z**2 / (beta**2 * A_brain)) *            (np.log(2 * me * c**2 * beta**2 * gamma**2 * Tmax / I_brain) - beta**2)
    return dEdx*rho_brain                                         

def simulate_proton_beam(E, num_protons=2000, target_size=20, beam_offset=10): 
    depths_energy_deposition = []
    E_threshold = 0.001  # Energy threshold in MeV below which the simulation stops

    
    for _ in range(num_protons):
        depth = beam_offset
        E_current = E
        while E_current > E_threshold and depth < (beam_offset + target_size):
            dEdx = bethe_bloch_energy(E_current)
            mean_step_length = 1 / dEdx
            step_length = -np.log(np.random.uniform()) * mean_step_length
            mean_delta_E = dEdx * step_length
            delta_E = np.random.normal(mean_delta_E, .05 * mean_delta_E)

            depth += step_length
            if beam_offset <= depth <= (beam_offset + target_size):
                # Convert delta_E from MeV/cm to MeV/cm^3 by dividing by the cross-sectional area
                energy_density = delta_E / cross_sectional_area
                depths_energy_deposition.append((depth - beam_offset, energy_density))

            E_current -= delta_E
            if E_current < 0:
                E_current = 0

    return depths_energy_deposition

#NEW
def read_csv_data(filename):
    depths = []
    energy_depositions = []
    with open(filename, 'r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip the header row
        for row in csv_reader:
            depths.append(float(row[0]))
            energy_depositions.append(float(row[1]))
    return depths, energy_depositions
#END

# Initialize plot
plt.figure(figsize=(12, 7))

#NEW
energies = np.arange(50, 251, 50)
for i, E in enumerate(energies):
    filename = f'proton_beam_{E}_MeV_energy_deposition.csv'
    depths, energy_deposits = read_csv_data(filename)
    plt.plot(depths, energy_deposits, label=f'W {E} MeV', color=water_colors[i], alpha=0.3)
#END

# Simulate and plot for each energy level
energies = np.arange(50, 251, 50)  
for i, E in enumerate(energies):
    result = simulate_proton_beam(E)
    depths, energy_deposits = zip(*result)
    depth_bins = np.linspace(0, max(depths), 100)
    dose_profile, _ = np.histogram(depths, bins=depth_bins, weights=energy_deposits)
    plt.plot(depth_bins[:-1], dose_profile, label=f'Br {E} MeV', color=colors[i])  # Use color from list

#edit the graph to look bette
plt.xlabel('Penetration Depth (cm)', fontsize=13, labelpad=10)
plt.ylabel('Energy Deposition (MeV/cm³)', fontsize=13, labelpad=10)
plt.title('Comparative Proton Beam Energy Deposition in Brain Tissue and Water', fontsize=14, pad=20)
plt.legend(title="Proton Energy", title_fontsize='large', ncol=2, fontsize='medium')
plt.tick_params(axis='both', which='major', labelsize=12)

# Use scientific notation for y-axis
ax = plt.gca()  # Get current axis
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
y_formatter = ScalarFormatter(useMathText=True) 
y_formatter.set_scientific(True) 
y_formatter.set_powerlimits((-3,3))
ax.yaxis.set_major_formatter(y_formatter)

plt.grid(True)
plt.tight_layout(pad=2.0)  # Adjust the padding between and around subplots
plt.show()


# In[70]:


#Molar Mass of Brain Tissue A_Brain

atomic_masses = {
    1: 1.008,   # Hydrogen
    6: 12.011,  # Carbon
    7: 14.007,  # Nitrogen
    8: 15.999,  # Oxygen
    11: 22.990, # Sodium
    12: 24.305, # Magnesium
    15: 30.974, # Phosphorus
    16: 32.065, # Sulfur
    17: 35.453, # Chlorine
    19: 39.098, # Potassium
    20: 40.078, # Calcium
    26: 55.845, # Iron
    30: 65.38   # Zinc
}

weight_fractions = {
    1: 0.110667,
    6: 0.125420,
    7: 0.013280,
    8: 0.737723,
    11: 0.001840,
    12: 0.000150,
    15: 0.003540,
    16: 0.001770,
    17: 0.002360,
    19: 0.003100,
    20: 0.000090,
    26: 0.000050,
    30: 0.000010
}

molar_mass_brain = sum(weight_fractions[el] * atomic_masses[el] for el in weight_fractions)
print(f"The estimated molar mass of brain tissue is {molar_mass_brain:.2f} g/mol")


# In[71]:


#Eff Atomic Number of Brain Tissue Z_brain


atomic_numbers = {
    1: 1,     # Hydrogen
    6: 6,     # Carbon
    7: 7,     # Nitrogen
    8: 8,     # Oxygen
    11: 11,   # Sodium
    12: 12,   # Magnesium
    15: 15,   # Phosphorus
    16: 16,   # Sulfur
    17: 17,   # Chlorine
    19: 19,   # Potassium
    20: 20,   # Calcium
    26: 26,   # Iron
    30: 30    # Zinc
}

weight_fractions = {
    1: 0.110667,
    6: 0.125420,
    7: 0.013280,
    8: 0.737723,
    11: 0.001840,
    12: 0.000150,
    15: 0.003540,
    16: 0.001770,
    17: 0.002360,
    19: 0.003100,
    20: 0.000090,
    26: 0.000050,
    30: 0.000010
}

p = 3

Z_eff = (sum(weight_fractions[Z] * (atomic_numbers[Z] ** p) for Z in atomic_numbers)) ** (1/p)

print(f"The effective atomic number for brain tissue is: {Z_eff:.2f}")


# In[22]:


#Water Vapor (https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=277)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

colors = [ '#000080', '#0000cc', '#0033ff', '#1e90ff', '#4169e1']

# Constants GOOD
me = 0.511                              # MeV/c^2, electron mass
c = 3e10                                # cm/s, speed of light
K = 0.307075                            # MeV cm^2 / mol, constant in Bethe-Bloch formula
m_proton = 938.27                       # MeV/c^2, proton mass
z = 1                                   # charge of proton

#Water Vapor Constants 
Z_water = 7.22                          # Effective atomic number for water (The physics of proton therapy Wayne D Newhauser and Rui Zhang))
A_water = 14.32                         # g/mol, molar mass of water
I_water = 71.6e-6                       # MeV, mean excitation energy for water
rho_water = 7.56182e-4                 # g/cm³, density of water

#Beam Constants
beam_radius = 0.5                                # radius of the proton beam in cm
cross_sectional_area = np.pi * beam_radius**2    # area in cm^2

#Math methods 

# Bethe-Bloch energy loss function
def bethe_bloch_energy(E_kin):
    beta = np.sqrt(1 - (m_proton / (E_kin + m_proton))**2)
    gamma = 1 / np.sqrt(1 - beta**2)
    Tmax = (2 * me * c**2 * beta**2 * gamma**2) /            (1 + (2 * gamma * me / m_proton) + (me / m_proton)**2)
    dEdx = K * Z_water * (z**2 / (beta**2 * A_water)) *            (np.log(2 * me * c**2 * beta**2 * gamma**2 * Tmax / I_water) - beta**2)
    return dEdx*rho_water                                         

def simulate_proton_beam(E, num_protons=500, target_size=20, beam_offset=10): 
    depths_energy_deposition = []
    E_threshold = 0.001  # Energy threshold in MeV below which the simulation stops

    
    for _ in range(num_protons):
        depth = beam_offset
        E_current = E
        while E_current > E_threshold and depth < (beam_offset + target_size):
            dEdx = bethe_bloch_energy(E_current)
            mean_step_length = 1 / dEdx
            step_length = -np.log(np.random.uniform()) * mean_step_length
            mean_delta_E = dEdx * step_length
            delta_E = np.random.normal(mean_delta_E, .001 * mean_delta_E)

            depth += step_length
            if beam_offset <= depth <= (beam_offset + target_size):
                # Convert delta_E from MeV/cm to MeV/cm^3 by dividing by the cross-sectional area
                energy_density = delta_E / cross_sectional_area
                depths_energy_deposition.append((depth - beam_offset, energy_density))

            E_current -= delta_E
            if E_current < 0:
                E_current = 0

    return depths_energy_deposition

# Initialize plot
plt.figure(figsize=(10, 6))


# Simulate and plot for each energy level
energies = np.arange(50, 251, 50)  
for i, E in enumerate(energies):
    result = simulate_proton_beam(E)
    depths, energy_deposits = zip(*result)
    depth_bins = np.linspace(0, max(depths), 100)
    dose_profile, _ = np.histogram(depths, bins=depth_bins, weights=energy_deposits)
    plt.plot(depth_bins[:-1], dose_profile, label=f'{E} MeV', color=colors[i])  # Use color from list


#edit the graph to look bette
plt.xlabel('Penetration Depth (cm)', fontsize=13, labelpad=10)
plt.ylabel('Energy Deposition (MeV/cm³)', fontsize=13, labelpad=10)
plt.title('Proton Beam Energy Deposition in Water Vapor', fontsize=14, pad=20)
plt.legend(title="Proton Energy", fontsize='large', title_fontsize='large')
plt.tick_params(axis='both', which='major', labelsize=12)

# Use scientific notation for y-axis
ax = plt.gca()                                                 # Get current axis
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
y_formatter = ScalarFormatter(useMathText=True) 
y_formatter.set_scientific(True) 
y_formatter.set_powerlimits((-3,3))
ax.yaxis.set_major_formatter(y_formatter)

plt.grid(True)
plt.tight_layout(pad=2.0)                                      # Adjust the padding between and around subplots
plt.show()


# In[ ]:


#Check https://www.sciencedirect.com/science/article/abs/pii/0092640X84900020?via%3Dihub for density correction

