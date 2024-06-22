import os
import pymser
import pandas as pd

from modules.raspa_input import create_GCMC_input
from modules.raspa_output import parse_GCMC

# Ignore warnings
import warnings
warnings.filterwarnings("ignore")

# Define the conditions of the simulation

OutputPath = 'GCMC'
FrameworkName = 'MgMOF-74'
UnitCells = 'auto'
NumberOfCycles = 1000
UseChargesFromCIFFile = 'yes'
GasComposition = {'CO2': 1.0}
ExternalTemperature = 298.0
ExternalPressure = 1e4
DesiredProdCycles = 2000

# --------------------------------------------

# Create a custom header for print statements

header = f"""
    ============================================================
             Automatic GCMC Simulation with pyMSER
    ============================================================
    Framework Name        : {FrameworkName}
    External Temperature  : {ExternalTemperature} K
    External Pressure     : {ExternalPressure} Pa
    Gas Composition       : {GasComposition}
    Desired Prod. Cycles  : {DesiredProdCycles}
    Output Path           : {OutputPath}
    ===========================================================

    Running RASPA simulation...
    """

print(header)

os.chdir(OutputPath)

equilibrated = False
nstep = 0
maxSteps = 20

while not equilibrated and nstep < maxSteps:

    nstep += 1
    print(f'       > Running iteration {nstep}...')

    create_GCMC_input(
        path='.',
        FrameworkName=FrameworkName,
        UnitCells=UnitCells,
        NumberOfCycles=NumberOfCycles,
        ForceField='local',
        UseChargesFromCIFFile=UseChargesFromCIFFile,
        GasComposition=GasComposition,
        ExternalTemperature=ExternalTemperature,
        ExternalPressure=ExternalPressure,
        RestartFile='no' if nstep == 1 else 'yes',
        )

    os.system('${RASPA_DIR}/bin/simulate simulation.input > raspalog.txt 2>&1')

    parse_GCMC(
        output_folder='Output/System_0',
        FrameworkName=FrameworkName,
        GasComposition=GasComposition,
        ExternalTemperature=ExternalTemperature,
        ExternalPressure=ExternalPressure,
        NumberOfCycles=NumberOfCycles,
        PrintEvery=1
        )

    dataFrame = pd.read_csv(f'Output/System_0/raspa_{ExternalTemperature:.6f}_{ExternalPressure}.csv')

    eqDict = pymser.equilibrate(dataFrame['N_ads'], print_results=False)

    equilibrated = len(dataFrame['N_ads']) - eqDict['t0'] > DesiredProdCycles

    if equilibrated:

        log_text = '==============================================================================\n'

        print("       > Success! Found {} production cycles. Analyzing final data...\n\n".format(
            len(dataFrame['N_ads']) - eqDict['t0']))

        eqDict = pymser.equilibrate(dataFrame['N_ads'], print_results=True)

        log_text = '==============================================================================\n'

        for i, gas in enumerate(GasComposition.keys()):
            eq_data = pymser.calc_equilibrated_average(
                data=dataFrame[f'{gas}_[mol/kg]'],
                eq_index=eqDict['t0'],
                uncertainty='uSD',
                ac_time=eqDict['ac_time']
                )

            enthalpy_data = pymser.calc_equilibrated_enthalpy(
                energy=dataFrame['total_[K]'],
                number_of_molecules=dataFrame[f'{gas}_[molecules/uc]'],
                temperature=ExternalTemperature,
                eq_index=eqDict['t0'],
                uncertainty='uSD',
                ac_time=eqDict['ac_time'])

            log_text += f'Component {i} [{gas}]\n'
            log_text += '-'*60 + '\n'
            log_text += 'Average loading absolute [mol/kg framework] {:20.10f} +/-  {:20.10f}\n'.format(*eq_data)
            log_text += 'Total Enthalpy of adsorption [KJ/MOL]       {:20.10f} +/-  {:20.10f}\n'.format(*enthalpy_data)
            log_text += '==============================================================================\n'

            print(log_text)
            with open(f'{FrameworkName}_{ExternalTemperature:.6f}_{ExternalPressure}.log', 'a') as f:
                f.write(log_text)
    else:
        print("       > Found only {}/{} production cycles. Running more {} cycles.".format(
            len(dataFrame['N_ads']) - eqDict['t0'], DesiredProdCycles, NumberOfCycles)
            )

        os.makedirs('RestartInitial/System_0', exist_ok=True)

        # Copy the file from Restart/System_0 to RestartInitial/System_0
        os.system('cp -r Restart/System_0/* RestartInitial/System_0/')
