import gemmi
from glob import glob
import re
import numpy as np

def calculate_Perpendicular_Widths(cif_filename: str) -> tuple:
    """
    Calculate the perpendicular widths of the unit cell.
    RASPA considers the perpendicular directions as the directions perpendicular to the `ab`,
    `bc`, and `ca` planes. Thus, the directions depend on the crystallographic vectors `a`, `b`,
    and `c`.
    The length in the perpendicular directions are the projections of the crystallographic vectors
    on the vectors `a x b`, `b x c`, and `c x a`. (here `x` means cross product)
    """
    # Read data from CIF file
    cif = gemmi.cif.read_file(cif_filename).sole_block()
    a = float(cif.find_value('_cell_length_a').split('(')[0])
    b = float(cif.find_value('_cell_length_b').split('(')[0])
    c = float(cif.find_value('_cell_length_c').split('(')[0])
    beta = float(cif.find_value('_cell_angle_beta').split('(')[0]) * np.pi / 180.0
    gamma = float(cif.find_value('_cell_angle_gamma').split('(')[0]) * np.pi / 180.0
    alpha = float(cif.find_value('_cell_angle_alpha').split('(')[0]) * np.pi / 180.0

    # Calculate the nu value
    nu = (np.cos(alpha) - np.cos(gamma) * np.cos(beta)) / np.sin(gamma)

    # Build the transformation matrix as a numpy array
    CellBox = np.array([[a, 0.0, 0.0],
                        [b * np.cos(gamma), b * np.sin(gamma), 0.0],
                        [c * np.cos(beta), c * nu, c * np.sqrt(1.0 - np.cos(beta)**2 - nu**2)]])

    # Calculate the cross products
    axb = np.cross(CellBox[0], CellBox[1])
    bxc = np.cross(CellBox[1], CellBox[2])
    cxa = np.cross(CellBox[2], CellBox[0])

    # Calculates the volume of the unit cell
    V = np.dot(np.cross(CellBox[0], CellBox[1]), CellBox[2])

    # Calculate perpendicular widths
    p_width_1 = V / np.linalg.norm(bxc)
    p_width_2 = V / np.linalg.norm(cxa)
    p_width_3 = V / np.linalg.norm(axb)

    return p_width_1, p_width_2, p_width_3


def calculate_UnitCells(cif_filename: str, cutoff: float) -> str:
    """
    Calculate the number of unit cell repetitions so that all supercell lengths are larger than
    twice the interaction potential cut-off radius.
    """

    # Calculate the perpendicular widths
    p_width_1, p_width_2, p_width_3 = calculate_Perpendicular_Widths(cif_filename)

    # Calculate UnitCells string
    uc_array = np.ceil(2.0 * cutoff / np.array([p_width_1, p_width_2, p_width_3])).astype(int)
    unit_cells = ' '.join(map(str, uc_array))

    return unit_cells


def get_pseudoatoms(molecule: str) -> list:
    """
    Returns the pseudoatoms of a given molecule.
    If the molecule is not in the supported list will returns `None`.
    Parameters
    ----------
    molecule : string
        Molecule name. Could be CO2, N2, O2, or H2O.
    Returns
    ----------
    pseudoatoms : list
        List containing the strings with the pseudotoms.
    """

    pseudoatoms_dict = {'CO2': ['C_co2', 'O_co2'],
                        'N2': ['N_n2', 'N_com'],
                        'O2': ['O_o2', 'O_com'],
                        'H2': ['H_h2', 'H_com'],
                        'CH4': ['CH4'],
                        'CO': ['C_co', 'CO_com', 'O_co'],
                        'H2O': ['Ow', 'Hw', 'Lw']}

    if molecule in list(pseudoatoms_dict.keys()):
        return pseudoatoms_dict[molecule]
    else:
        return None


def get_conversion_factors(output_folder: str,
                           FrameworkName: str,
                           ExternalTemperature: float,
                           ExternalPressure: float,):
    """
    Get the conversion factors for the units in the RASPA simulation.

    Parameters
    ----------
    path : string
        Path to the folder containing the RASPA output file.
    filename : string
        Name of the RASPA output file.
    Returns
    ----------
    conversion_factors : dict
        Dictionary containing the conversion factors.
    """

    # Read file into string
    filename = glob('{0}/output_{1}_*_{2:.6f}_{3:g}.data'.format(output_folder,
                                                                 FrameworkName,
                                                                 ExternalTemperature,
                                                                 ExternalPressure))[0]

    pattern = re.compile(r'Conversion factor molecules/unit cell -> (.+?):\s+(\d+\.\d+)')

    with open(filename, 'r') as f:
        lines = f.readlines()

    conversion_factors = {
        'mol/kg': [],
        'mg/g': [],
        'cm^3 STP/gr': [],
        'cm^3 STP/cm^3': []
    }
    for line in lines:
        if 'Conversion factor molecules/unit cell' in line:
            match = re.search(pattern, line)

            conversion_factors[match.group(1)].append(float(match.group(2)))
    
    return conversion_factors
