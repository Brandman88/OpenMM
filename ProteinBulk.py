import os
#import sys
#print(sys.path)
import hoomd
import hoomd.md
import gsd.hoomd
import numpy as np
import pandas as pd

#parameters='multi_run.dat'
#num_run is first number, says how many runs are we doing?
#start is the second number, says where we are in that list
#marker is are we at the start or end of a task. default is S for starting

#dest_file='1LC.txt'

cur_dir=os.getcwd()
dir_comp=f'{cur_dir}/completed_run'

def find_Protein_length(dest_file_OLC='OneLetterCode.dat'):
    seq_length = len(np.loadtxt('ThreeLetterCode.dat', dtype="U", unpack=True))
    return seq_length
        
def read_entire_csv_return_dict(dest_required_file_csv = "data.csv"):
    '''Define function to read entire CSV and return a dictionary.'''
    # Using pandas to read the CSV file located at the destination provided.
    df = pd.read_csv(dest_required_file_csv)

    # Converting the dataframe into a list of dictionaries where each dictionary represents a row of data.
    data = df.to_dict('records')

    # Return the list of dictionaries
    return data

def grab_info_from_csv_dict_row(start):
    '''Define function to grab information from a specific row in the CSV.'''
    # Call the function to read the CSV file and store the data in csv_dictionary
  
    
    csv_dictionary=read_entire_csv_return_dict()

    # Adjusting the row number to align with zero-indexed lists.
    csv_row_num=start-1 

    # Getting the information from the specific row in the CSV.
    current_row_info=csv_dictionary[csv_row_num]

    # Return the row information
    return current_row_info

def row_specific_info_csv(column_name):
    '''Define function to grab specific column information from a specific row in the CSV.'''
    num_run,start,marker,clean_list=read_multi_run()
    
    # Call the function to get the information of the specific row
    current_row_info=grab_info_from_csv_dict_row(start)

    # Get the information of the specific column from the row.
    specific_csv_row_info=current_row_info[column_name]

    # Return the specific column information
    return specific_csv_row_info

# This function reads a multi-run parameter file and extracts the parameters from it.
def read_multi_run(parameters='multi_run.dat'):
    # Open the specified parameter file.
    with open(parameters,'r') as para:
        lines = []
        # Read the file line by line and store each line in a list.
        for line in para:
            lines.append(line)
    para.close()
    num_run_raw=lines[0]
    start_raw=lines[1]
    marker_raw=lines[2]
    
    # Extract the number of proteins to run.
    if '#Number of proteins you care to run' in num_run_raw:
        num_run=int(num_run_raw.replace('#Number of proteins you care to run','').strip())
    else:
        num_run=int(num_run_raw)

    # Extract the start value.
    if '#Start value (if you are just starting the process should be on 1 by default)' in start_raw:
        start=int(start_raw.replace('#Start value (if you are just starting the process should be on 1 by default)','').strip())
    else:
        start=int(start_raw)

    # Extract the marker.
    if '#Put S or E based on starting or end process (don\'t really need to touch), if its I then set the other 2 numbers and run the file again.' in marker_raw:
        marker=marker_raw.replace('#Put S or E based on starting or end process (don\'t really need to touch), if its I then set the other 2 numbers and run the file again.','')
    else:
        marker=marker_raw
    marker=marker.strip()

    # Remove the lines that have been processed.
    lines.remove(num_run_raw)
    lines.remove(start_raw)
    lines.remove(marker_raw)
  
    clean_list=[]
  
    # Go through the remaining lines and clean them up.
    for line in lines:
        line=line.strip()
        # If the line is not empty,
        if line.strip():
            # If the line contains a newline character, remove it.
                if '\n' in line:
                    # Remove newline character and any leading or trailing white spaces
                    line=line.replace('\n','').strip()
                    clean_list.append(line)
                else:
                    # Add the line into the clean list after removing leading or trailing white spaces
                    clean_list.append(line.strip())
                
    # Return the number of runs, start value, marker, and the cleaned list of lines
    return num_run,start,marker,clean_list

def box_size_equation():
    est_box_size = int(((4.08744 * 10**-3) * find_Protein_length()**2.35278) + 207)
    return est_box_size

box_size=box_size_equation()
hoomd.context.initialize("");

# ========================= System Parameters =======================================
ParticleN=find_Protein_length() #imports number of ParticleN
dimen = 3
# Box Dimentions
BOX_L = box_size

#kappa = 0.0

# ========================= Simulation Constants =======================================
T_Kelvin=float(row_specific_info_csv('Absolute Temperature'))
# In kJ units
#KT = T_Kelvin*8.3145e-3 # Boltzmann constant in kJ/mol/K
#EPSILON = 0.2 * 4.184 # kJ/mol

# In kCal units
KT = T_Kelvin*0.001987204259 # Boltzmann constant in kCal/mol/K
EPSILON = 0.2 # KCal/mol

ionic_concentration = row_specific_info_csv('Ionic Concentration')

if isinstance(ionic_concentration, int):
    ionic_concentration = int(ionic_concentration)  # Convert to int if it's already an integer
elif isinstance(ionic_concentration, float):
    ionic_concentration = float(ionic_concentration)  # Convert to float if it's already a float
else:
    try:
        ionic_concentration = float(eval(ionic_concentration))  # Evaluate the string expression and convert to float
    except (SyntaxError, NameError, TypeError, ValueError):
        # Handle the case where the value cannot be evaluated or converted
        # Set a default value or raise an exception, depending on your requirement
        ionic_concentration = None  # or set a default value

# Now you have the ionic_concentration variable as an integer or float, or None if it couldn't be evaluated or converted.

fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3 #temperature dependent dielectric constant of water	
epsw = fepsw(T_Kelvin) # dielectric constant of water at T 	
lB = (1.6021766**2/(4*np.pi*8.854188*epsw))*(6.022*1000/KT)/4.184 # Bjerrum length in nm	
#lB = (1.6021766**2/(4*np.pi*8.854188*80))*(6.022*1000/KT)/4.184 #  Bjerrum length in nm	

yukawa_eps = lB*KT	
yukawa_kappa = np.sqrt(8*np.pi*lB*ionic_concentration*6.022/10)	

# ========================= System Parameters =======================================
fi_sys = open("System_Parameters.dat", "w")
fi_sys.write("Parameters and Constants ============================\n")
fi_sys.write("Number of particles: %d\n" % ParticleN)
fi_sys.write("Box length (nm): %f\n" % BOX_L)
fi_sys.write("Temperature (K): %d\n" % T_Kelvin)
fi_sys.write("LJ epsilon (kCal/mol): %f\n" % EPSILON)
fi_sys.write("Gas constant*T (kCal/mol): %f\n" % KT)
fi_sys.write(f"Dielectric constant of water at T={T_Kelvin}: %f\n" % epsw)
fi_sys.write("Bjerrum length (nm): %f\n" % lB)
fi_sys.write("Ionic concentration (M): %f\n" % ionic_concentration)
fi_sys.write("Yukawa epsilon (nm kCal/mol): %f\n" % yukawa_eps)
fi_sys.write("Yukawa kappa (nm^-1): %f\n" % yukawa_kappa)
fi_sys.write("=================================================== \n")
fi_sys.close()
# ====================================================================================

TIME_STEP = 0.01 # in picoseconds given current units
nwrite = 1000 # Frequency of saving configurations

Rouse_Steps = int((5*(ParticleN)**2.2)/TIME_STEP)
Run_Steps = Rouse_Steps
#Rouse_Steps = 1000000
#Run_Steps = 20000000

seed = np.random.randint(5000)
print ("Seed for the Run:", seed)

# ========================= Particles Connection Initialization =======================================
bond_gr=[]
for i in range(ParticleN-1):
	a_bond=[]
	a_bond.append(i)
	a_bond.append(i+1)
	bond_gr.append(a_bond)

angle_gr=[]	
for i in range(ParticleN-2):
	a_angle=[]
	a_angle.append(i)
	a_angle.append(i+1)
	a_angle.append(i+2)
	angle_gr.append(a_angle)
	
# ========================= System Initialization =======================================
system=gsd.hoomd.Snapshot()

aakeys = np.loadtxt('stats_module.dat', dtype=str, usecols=(0), unpack=True)
aakeys_mass = np.loadtxt('stats_module.dat', usecols=(1), unpack=True)
aakeys_chgs = np.loadtxt('stats_module.dat', usecols=(2), unpack=True)
aakeys_sigmas = np.loadtxt('stats_module.dat', usecols=(3), unpack=True)
aakeys_sigmas_scaled = aakeys_sigmas/10. # Angstrom to nm
aakeys_lambdas = np.loadtxt('stats_module.dat', usecols=(4), unpack=True)
aavalues = np.loadtxt('chain_param.dat', usecols=(1,2,3,4,5), unpack=True)

# Box description
system.configuration.dimensions = dimen
system.configuration.box = [BOX_L,BOX_L,BOX_L,0,0,0]

# Particle description
system.particles.N = ParticleN
system.particles.position = np.loadtxt(str(ParticleN))
system.particles.types = aakeys
system.particles.typeid = aavalues[0].astype(int)
system.particles.mass = aavalues[1]
system.particles.charge = aavalues[2]

# Bond description
system.bonds.N = ParticleN-1
system.bonds.types = ['AA_bond']
system.bonds.group = bond_gr

# Angle description
system.angles.N = ParticleN-2
system.angles.types = ['AA_angle']
system.angles.group = angle_gr

# Write intial chain config into gsd file
fi_start = gsd.hoomd.open(name='Start_Config.gsd', mode='wb')
fi_start.append(system)
#fi_start.close()

system = hoomd.init.read_gsd('Start_Config.gsd')

# ========================= Neighbor List =======================================
#nl = hoomd.md.nlist.cell(); # Change to tree for large box size
nl = hoomd.md.nlist.tree();

# ========================= Force Fields =======================================
# Custom Ashbaugh Potential
lj1 = hoomd.md.pair.lj(r_cut=3.5, nlist=nl, name="1")
lj2 = hoomd.md.pair.lj(r_cut=3.5, nlist=nl, name="2")

# Harmonic Bonds
harmonic=hoomd.md.bond.harmonic()
#harmonic.bond_coeff.set('AA_bond', k=8033, r0=0.38) # k in kJ/mol/nm^2 and r0 in nm
harmonic.bond_coeff.set('AA_bond', k=1920, r0=0.38) # k in kCal/mol/nm^2 and r0 in nm


# FENE Potential
#fene = hoomd.md.bond.fene()
#fene.bond_coeff.set('AA_bond', k=30.0, r0=1.5, sigma=1.0, epsilon= 0.0)
		
# Angle Potential: 
# V = k[1-cos(theta - theta0)], where theta0 = np.pi
# Force: T = - dV/d(theta)

#def bend_pot(theta, kappa):
#	V = kappa * (1.0+np.cos(theta));
#	T = kappa*np.sin(theta);
#	return (V,T)

#btable = hoomd.md.angle.table(width=1000)
#btable.angle_coeff.set('AA_angle', func=bend_pot, coeff=dict(kappa=kappa))

# Electrostatics
yukawa = hoomd.md.pair.yukawa(r_cut=3.5, nlist=nl)

# ========================= Pair Coefficients =======================================
for i in range (len(aakeys)):
	for j in range (len(aakeys)):
		sig_eff = 0.5*(aakeys_sigmas_scaled[i] + aakeys_sigmas_scaled[j]) 
		lambda_eff = 0.5*(aakeys_lambdas[i] + aakeys_lambdas[j]) 
		#print (i, j, aakeys[i], aakeys[j], sig_eff, lambda_eff)

		lj1.pair_coeff.set(aakeys[i], aakeys[j], epsilon=EPSILON*(1-lambda_eff), sigma=sig_eff, r_cut=(2**(1./6))*sig_eff)
		lj2.pair_coeff.set(aakeys[i], aakeys[j], epsilon=EPSILON*lambda_eff, sigma=sig_eff, r_cut=3.5)
  
		yukawa.pair_coeff.set(aakeys[i], aakeys[j], epsilon=yukawa_eps*aakeys_chgs[i]*aakeys_chgs[j], kappa=yukawa_kappa, r_cut=3.5)

lj1.set_params(mode='shift')
yukawa.set_params(mode='shift')
nl.reset_exclusions(exclusions = ['bond']);

# ========================= MD Integrator =======================================
hoomd.md.integrate.mode_standard(dt=TIME_STEP);
all = hoomd.group.all();
integrator = hoomd.md.integrate.langevin(group=all, kT=KT, seed=seed, noiseless_t=False, noiseless_r=False)

for i,aa in enumerate(aakeys):
	GAMMA = aakeys_mass[i]/1000.
	integrator.set_gamma_r(aa, gamma_r=GAMMA)
	
# ========================= Warmup Run =======================================
hoomd.run(Rouse_Steps, quiet=True);
hoomd.dump.gsd("Equilibration_Config.gsd", period=nwrite, group=all, overwrite=True);
integrator.disable()

# ========================= MD Integrator =======================================
hoomd.md.integrate.mode_standard(dt=TIME_STEP);
all = hoomd.group.all();
integrator = hoomd.md.integrate.langevin(group=all, kT=KT, seed=seed, noiseless_t=False, noiseless_r=False)
for i,aa in enumerate(aakeys):
	GAMMA = aakeys_mass[i]/1000.
	integrator.set_gamma_r(aa, gamma_r=GAMMA)

# ========================= Print Thrmodynamic Quantities =======================================
#hoomd.analyze.log(filename="mddata.dat",quantities=['potential_energy','kinetic_energy','pair_lj_energy','bond_fene_energy', 'temperature'],period=10*nwrite,overwrite=True);
#hoomd.analyze.log(filename="mddata.dat",quantities=['potential_energy','kinetic_energy','pair_ashbaugh_energy','bond_harmonic_energy','pair_yukawa_energy', 'temperature'],period=10*nwrite,overwrite=True);
hoomd.analyze.log(filename="mddata.dat",quantities=['potential_energy','kinetic_energy','bond_harmonic_energy','pair_yukawa_energy', 'temperature'],period=10*nwrite,overwrite=True);

# ========================= Trajectory Print for Movie =======================================
hoomd.dump.gsd("Running_Config.gsd", period=nwrite, group=all, overwrite=True);
#hoomd.dump.dcd(filename="Running_Config.dcd", period=nwrite)

# ========================= Run Steps =======================================
hoomd.run(Run_Steps);

# ========================= Run Analysis =======================================
#os.system("python3 config_analysis.py")

