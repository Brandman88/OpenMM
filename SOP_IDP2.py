from simtk.openmm import app
import simtk.openmm as omm
from simtk import unit
import time
import sys
import argparse
import build
import force
import pandas as pd

def read_entire_csv_return_dict(dest_required_file_csv = "data_multi.csv"):
    '''Define function to read entire CSV and return a dictionary.'''
    # Using pandas to read the CSV file located at the destination provided.
    df = pd.read_csv(dest_required_file_csv)
    df['letters_count'] = df['Sequence'].apply(len)
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

def get_equilibrium_data_forfeiture(parameters='multi_run.dat', filename='data_multi.csv'):
    # Call the read_multi_run function to get the start value
    num_run,start,marker,clean_list=read_multi_run()
    
    # Read the CSV file into a DataFrame
    df = pd.read_csv(filename)
    
    # Get the row index based on the start value
    row_index = start - 1
    
    # Get the value for the "Equilibrium Data Forfeiture" header
    equilibrium_data_forfeiture = df.loc[row_index, 'Equilibrium Data Forfeiture']
    
    # Check if the value is empty or greater than or equal to 1
    if equilibrium_data_forfeiture == "" or float(equilibrium_data_forfeiture) >= 100:
        equilibrium_data_forfeiture = 70
    
    return equilibrium_data_forfeiture,num_run,start,marker,clean_list

equilibrium_data_forfeiture,num_run,start,marker,clean_list = get_equilibrium_data_forfeiture()  # The % of the data that will be sacrificed to claim equilibrium 


def parse_arguments_from_csv(start,filename='data_multi.csv'):
    TIME_STEP=0.01 #in Picoseconds
    df = pd.read_csv(filename)
    ParticleN=row_specific_info_csv('letters_count')
    row_index = int(start) - 1  # Assuming start is a string representing the row number
    arguments = df.iloc[row_index].to_dict()

    # Convert NaN values to None
    arguments = {k: v if pd.notna(v) else None for k, v in arguments.items()}

    frequency = arguments['Frequency']
    frequency_default = 1000  # Set your desired default value here

    cutoff = arguments['Cutoff']
    cutoff_default = 40.0  # Set your desired default value here
    
    trajectory = arguments['Trajectory']
    trajectory_default = 'Running_Config.pdb'  # Set your desired default value here
    
    output = arguments['Output']
    output_default = 'Running_Config.out'  # Set your desired default value here

    number_of_steps = arguments['Number of Steps']
    number_of_steps_default = int((3*(ParticleN)**2.2)/TIME_STEP)  # Set your desired default value here

    
    # Convert frequency to int or use the default value
    frequency = int(frequency) if frequency is not None else frequency_default
    arguments['Frequency'] = frequency

    cutoff = float(cutoff) if cutoff is not None else cutoff_default
    arguments['Cutoff'] = cutoff
    
    trajectory = trajectory if trajectory is not None else trajectory_default
    arguments['Trajectory'] = trajectory
    
    output = output if output is not None else output_default
    arguments['Output'] = output
    
    number_of_steps = number_of_steps if number_of_steps is not None else number_of_steps_default
    arguments['Number of Steps'] = number_of_steps    
    
    return arguments



KELVIN_TO_KT = unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB / unit.kilocalorie_per_mole


filename = 'data_multi.csv'  # Replace with your CSV file name

arguments = parse_arguments_from_csv(start,filename)

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Coarse-grained SOP_IDP simulation using OpenMM')

# Add arguments from the CSV file
parser.add_argument('-f', '--sequence', type=str, help='input sequence', default=arguments['Sequence'])
parser.add_argument('-K', '--monovalent_concentration', type=float, default=float(arguments['Monovalent Concentration']),
                    help='Monovalent concentration (mM) [150.0]')
parser.add_argument('-c', '--cutoff', type=float, default=float(arguments['Cutoff']),
                    help='Cutoff distance for electrostatics (A) [40.0]')
parser.add_argument('-T', '--temperature', type=float, default=float(arguments['Absolute Temperature']),
                    help='Temperature (K) [293.15]')
parser.add_argument('-t', '--traj', type=str, default=arguments['Trajectory'],
                    help='trajectory output')
parser.add_argument('-o', '--output', type=str, default=arguments['Output'],
                    help='status and energy output')
parser.add_argument('-x', '--frequency', type=int, default=int(arguments['Frequency']),
                    help='output frequency')
parser.add_argument('-s', '--step', type=int, default=int(arguments['Number of Steps']),
                    help='Number of step [10000]')

args = parser.parse_args()

class simu:    ### structure to group all simulation parameter
    temp = 0.
    Kconc = 0.
    Nstep = 0
    epsilon = 0.
    cutoff = 40.

simu.temp = (args.temperature)*unit.kelvin
simu.Nstep = args.step
simu.Kconc = args.monovalent_concentration
simu.cutoff = args.cutoff

equilibrium_steps=simu.Nstep*(equilibrium_data_forfeiture/100)
steps_left=simu.Nstep-equilibrium_steps

T_unitless = simu.temp * KELVIN_TO_KT
print("T_unitless  ", T_unitless)
simu.epsilon = 296.0736276 - 619.2813716 * T_unitless + 531.2826741 * T_unitless**2 - 180.0369914 * T_unitless**3;
print("epsilon  ", simu.epsilon)
simu.l_Bjerrum = 332.0637*unit.angstroms / simu.epsilon
print("Bjerrum length  ", simu.l_Bjerrum / T_unitless)
simu.kappa = unit.sqrt (4*3.14159 * simu.l_Bjerrum * 2*simu.Kconc*6.022e-7 / (T_unitless * unit.angstrom**3))
print("kappa   ", simu.kappa)

forcefield = app.ForceField('SOP_IDP2.xml')
topology = None
positions = None

if args.sequence != None:
    print("Building from sequence %s ..." % args.sequence)
    topology, positions = build.build_by_seq(args.sequence, forcefield)
else:
    print("Need sequence !!!")
    sys.exit()

system = forcefield.createSystem(topology)

force.add_bond_force (topology, system, 0)
force.add_exV_force  (topology, system, 1)
force.add_DH_force   (topology, system, simu, 2)
force.add_statistical_force (topology, system, 3)
totalforcegroup = 3

integrator = omm.LangevinIntegrator(simu.temp, 0.5/unit.picosecond, 10*unit.femtoseconds)
simulation = app.Simulation(topology, system, integrator)

simulation.context.setPositions(positions)
print("Initial energy   %f   kcal/mol" % (simulation.context.getState(getEnergy=True).getPotentialEnergy() / unit.kilocalorie_per_mole))

state = simulation.context.getState(getPositions=True)
app.PDBFile.writeFile(topology, state.getPositions(), open("input.pdb", "w"), keepIds=True)

simulation.context.setVelocitiesToTemperature(simu.temp)
simulation.step(equilibrium_steps)
simulation.reporters.append(app.PDBReporter(args.traj, args.frequency))
simulation.reporters.append(app.StateDataReporter(args.output, args.frequency, elapsedTime=True, step=True, speed=True, progress=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, remainingTime=True, totalSteps=simu.Nstep, separator=','))

print('Running simulation ...')
t0 = time.time()
simulation.step(steps_left)
prodtime = time.time() - t0
print("Simulation speed: % .2e steps/day" % (86400*steps_left/(prodtime)))
