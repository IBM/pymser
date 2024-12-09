import os
import argparse
import json
import pymser
import numpy as np
import pandas as pd

from modules.raspa_input import create_GCMC_input
from modules.raspa_output import parse_GCMC
from modules.tools import get_conversion_factors

# Ignore warnings
import warnings
warnings.filterwarnings("ignore")


# Required parameters
parser = argparse.ArgumentParser(
    description='on-the-fly RASPA simulations with pyMSER')
parser.add_argument('output_folder',
                    type=str,
                    action='store',
                    metavar='OUTPUT_FOLDER',
                    help='Directory to save the files of the calculations')
parser.add_argument('--FrameworkName',
                    type=str,
                    required=True,
                    action='store',
                    metavar='FRAMEWORK_NAME',
                    help='Name of the framework to be simulated')
parser.add_argument('--ExternalPressure',
                    type=float,
                    required=True,
                    action='store',
                    metavar='PRESSURE_LIST',
                    help='External pressure in Pascal.')
# Optional parameters
parser.add_argument('--NumberOfProdCycles',
                    type=int,
                    required=False,
                    action='store',
                    metavar='NUMBER_OF_CYCLES',
                    default=5000,
                    help='Number of desired production cycles.')
parser.add_argument('--AddCycles',
                    type=int,
                    required=False,
                    action='store',
                    metavar='ADD_CYCLES',
                    default=1000,
                    help='Number of aditional cycles if NumberOfProdCycles was not achieved.')
parser.add_argument('--ExternalTemperature',
                    type=float,
                    required=False,
                    action='store',
                    metavar='TEMPERATURE',
                    default=298.0,
                    help='Temperature of the simulation in Kelvin')
parser.add_argument('--UnitCells',
                    type=str,
                    required=False,
                    action='store',
                    metavar='UNIT_CELLS',
                    default='auto',
                    help='Number of unit cells to simulate. Can be "auto" or "NX,NY,NZ" string.')
parser.add_argument('--UseChargesFromCIFFile',
                    default=False,
                    required=False,
                    action='store_true',
                    help='Use charges from CIF file.')
parser.add_argument('--GasComposition',
                    type=str,
                    required=False,
                    action='store',
                    metavar='GAS_COMPOSITION',
                    default='{"CO2": 1.0}',
                    help='Type of dispersion correction used for all calculations.')


arg = parser.parse_args()

arg.GasComposition = json.loads(arg.GasComposition)

header = f"""
==============================================================================
                    Automatic GCMC Simulation with pyMSER
==============================================================================
Framework Name        : {arg.FrameworkName}
External Temperature  : {arg.ExternalTemperature} K
External Pressure     : {arg.ExternalPressure} Pa
Gas Composition       : {arg.GasComposition}
Desired Prod. Cycles  : {arg.NumberOfProdCycles}
Output Path           : {arg.output_folder}
==============================================================================

    Running RASPA simulation...
    """

print(header)

os.chdir(arg.output_folder)

equilibrated = False
nstep = 0
maxSteps = 20

while not equilibrated and nstep < maxSteps:

    nstep += 1
    print(f'       > Running iteration {nstep}...')

    create_GCMC_input(
        path='.',
        FrameworkName=arg.FrameworkName,
        UnitCells=arg.UnitCells,
        NumberOfCycles=arg.AddCycles,
        ForceField='local',
        UseChargesFromCIFFile=arg.UseChargesFromCIFFile,
        GasComposition=arg.GasComposition,
        ExternalTemperature=arg.ExternalTemperature,
        ExternalPressure=arg.ExternalPressure,
        RestartFile='no' if nstep == 1 else 'yes',
        )

    os.system('${RASPA_DIR}/bin/simulate simulation.input > raspalog.txt 2>&1')

    parse_GCMC(
        output_folder='Output/System_0',
        FrameworkName=arg.FrameworkName,
        GasComposition=arg.GasComposition,
        ExternalTemperature=arg.ExternalTemperature,
        ExternalPressure=arg.ExternalPressure,
        NumberOfCycles=arg.AddCycles,
        PrintEvery=1
        )

    csv_file = f'Output/System_0/raspa_{arg.ExternalTemperature:.6f}_{arg.ExternalPressure}.csv'
    dataFrame = pd.read_csv(csv_file)

    eqDict = pymser.equilibrate(dataFrame['N_ads'], print_results=False)

    equilibrated = len(dataFrame['N_ads']) - eqDict['t0'] > arg.NumberOfProdCycles

    if equilibrated:

        log_text = '=============================================================================\n'

        print("       > Success! Found {} production cycles. Analyzing final data...\n\n".format(
            len(dataFrame['N_ads']) - eqDict['t0']))

        eqDict = pymser.equilibrate(dataFrame['N_ads'], print_results=True)

        log_text = '=============================================================================\n'

        convFactors = get_conversion_factors(
            output_folder='Output/System_0',
            FrameworkName=arg.FrameworkName,
            ExternalTemperature=arg.ExternalTemperature,
            ExternalPressure=arg.ExternalPressure)

        for i, gas in enumerate(arg.GasComposition.keys()):
            eq_data = pymser.calc_equilibrated_average(
                data=dataFrame[f'{gas}_[molecules/uc]'],
                eq_index=eqDict['t0'],
                uncertainty='uSD',
                ac_time=eqDict['ac_time']
                )

            eq_data = np.array(eq_data)

            enthalpy_data = pymser.calc_equilibrated_enthalpy(
                energy=dataFrame['total_[K]'],
                number_of_molecules=dataFrame[f'{gas}_[N_ads]'],
                temperature=arg.ExternalTemperature,
                eq_index=eqDict['t0'],
                uncertainty='uSD',
                ac_time=int(arg.NumberOfProdCycles/5))

            log_text += f'Component {i} [{gas}]\n'
            log_text += '-'*75 + '\n'
            log_text += 'Average loading absolute [molecules/unit cell] {:20.10f} +/-  {:20.10f}\n'\
                .format(*eq_data)
            log_text += 'Average loading absolute [mol/kg framework]    {:20.10f} +/-  {:20.10f}\n'\
                .format(*eq_data * convFactors['mol/kg'][i])
            log_text += 'Average loading absolute [mg/g framework]      {:20.10f} +/-  {:20.10f}\n'\
                .format(*eq_data * convFactors['mg/g'][i])
            log_text += 'Average loading absolute [cm^3 STP/gr]         {:20.10f} +/-  {:20.10f}\n'\
                .format(*eq_data * convFactors['cm^3 STP/gr'][i])
            log_text += 'Average loading absolute [cm^3 STP/cm^3]       {:20.10f} +/-  {:20.10f}\n'\
                .format(*eq_data * convFactors['cm^3 STP/cm^3'][i])
            log_text += 'Enthalpy of adsorption [KJ/mol]                {:20.10f} +/-  {:20.10f}\n'\
                .format(*enthalpy_data)
            log_text += '========================================================================\n'

        print(log_text)
        log = f'{arg.FrameworkName}_{arg.ExternalTemperature:.6f}_{arg.ExternalPressure}.log', 'a'
        with open(log) as f:
            f.write(log_text)
    else:
        print("       > Found only {}/{} production cycles. Running more {} cycles.".format(
            len(dataFrame['N_ads']) - eqDict['t0'], arg.NumberOfProdCycles, arg.AddCycles)
            )

        os.makedirs('RestartInitial/System_0', exist_ok=True)

        # Copy the file from Restart/System_0 to RestartInitial/System_0
        os.system('cp -r Restart/System_0/* RestartInitial/System_0/')
