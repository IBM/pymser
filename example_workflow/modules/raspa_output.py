import os
import math
from RASPA2 import parse
from glob import glob


def parse_GCMC(output_folder: str,
               FrameworkName: str,
               ExternalTemperature: float,
               ExternalPressure: float,
               GasComposition: dict,
               NumberOfCycles: int,
               PrintEvery: int):
    # Read file into string
    input_file_name = glob('{0}/output_{1}_*_{2:.6f}_{3:g}.data'.format(output_folder,
                                                                        FrameworkName,
                                                                        ExternalTemperature,
                                                                        ExternalPressure))[0]
    with open(os.path.join(input_file_name), 'r') as f:
        raspa_string = f.read()

    # Parse string into dictionary and list
    raspa_dict = parse(raspa_string)
    raspa_list = raspa_string.split('\n')

    # Extract number of unit cells in the supercell
    unit_cells = [int(line.split(':')[1]) for line in raspa_list if 'Number of unitcells' in line]

    # Calculate mol/kg conversion factor for the supercell
    to_mol_kg = raspa_dict['MoleculeDefinitions']['Conversion factor molecules/unit cell -> mol/kg'][0]
    to_mol_kg /= math.prod(unit_cells)

    if os.path.isfile(os.path.join(output_folder, f'raspa_{ExternalTemperature:.6f}_{ExternalPressure}.csv')):
        append_data = True
        # Open the file
        with open(os.path.join(output_folder, f'raspa_{ExternalTemperature:.6f}_{ExternalPressure}.csv'), 'r') as f:
            lines = f.readlines()
            # Check if the last line is the same as the last cycle
            last_line = lines[-1].split(',')
            base_cycle = int(last_line[0])
            base_step = int(last_line[1])

    else:
        append_data = False
        base_cycle = 0
        base_step = 0

    # Build header string
    if not append_data:
        header = 'cycle,step,N_ads'
        for component in range(len(GasComposition)):
            cycle_key = f'Current cycle: 0 out of {NumberOfCycles}'
            component_key = f'Component {component}'
            molecule_name = raspa_dict[cycle_key][component_key][0]
            header += (
                f',{molecule_name}_[N_ads]'
                f',{molecule_name}_[molecules/uc]'
                f',{molecule_name}_[mol/kg]'
            )

        header += (
            ',total_[K]'
            ',host-host_[K]'
            ',host-adsorbate_[K]'
            ',host-cation_[K]'
            ',adsorbate-adsorbate_[K]'
            ',cation-cation_[K]'
            ',adsorbate-cation_[K]'
        )

        csv_output = header + '\n'
    else:
        csv_output = ''

    # For each cycle
    steps = base_step
    for cycle in range(0, NumberOfCycles, PrintEvery):
        cycle_key = f'Current cycle: {cycle} out of {NumberOfCycles}'
        number_of_adsorbates = int(raspa_dict[cycle_key]['Number of Adsorbates'][0])
        steps += max(20, number_of_adsorbates)
        line = (
            f'{cycle + base_cycle},'
            f' {steps},'
            f' {number_of_adsorbates}'
        )

        # For each component
        for component in range(len(GasComposition)):
            component_key = f'Component {component}'
            number_of_molecules = int(raspa_dict[cycle_key][component_key][2].split('/')[0])
            line += (
                f', {number_of_molecules:7}'
                f', {number_of_molecules / math.prod(unit_cells):7}'
                f', {number_of_molecules * to_mol_kg:.7f}'
            )

        for energy_term in [
            'Current total potential energy',
            'Current Host-Host energy',
            'Current Host-Adsorbate energy',
            'Current Host-Cation energy',
            'Current Adsorbate-Adsorbate energy',
            'Current Cation-Cation energy',
            'Current Adsorbate-Cation energy'
        ]:

            line += f',{raspa_dict[cycle_key][energy_term][0]:.7f}'

        csv_output += line + '\n'

    # Write string into file
    output_file_name = f'raspa_{ExternalTemperature:.6f}_{ExternalPressure}.csv'
    with open(os.path.join(output_folder, output_file_name), 'a') as f:
        f.write(csv_output)
