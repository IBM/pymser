import os
from textwrap import dedent
from modules.tools import calculate_UnitCells, get_pseudoatoms


def create_GCMC_input(path: str, FrameworkName: str, **kwargs):
    """
    Create the RASPA GCMC simulation input.
    Parameters
    ----------
    path : string
        Path where the file will be saved.
    FrameworkName : string
        Name of the structure. Must be the same name in the `.cif` file.
    """

    # Calculation parameters dictionary
    CALC_DICT = {
        'FrameworkName': FrameworkName,                                 # string
        'NumberOfCycles': 10000,                                        # int
        'NumberOfInitializationCycles': 0,                              # int
        'PrintEvery': 1,                                                # int
        'PrintPropertiesEvery': 1,                                      # int
        'ForceField': 'local',                                          # string
        'CutOffVDW': 12.8,                                              # float
        'CutOffChargeCharge': 12.8,                                     # float
        'CutOffChargeBondDipole': 12.8,                                 # float
        'CutOffBondDipoleBondDipole': 12.8,                             # float
        'EwaldPrecision': 1.0e-6,                                       # float
        'HeliumVoidFraction': 0.0,                                      # float
        'ExternalTemperature': 298.15,                                  # int
        'ExternalPressure': '100000',                                   # float or csv
        'UseChargesFromCIFFile': 'yes',                                 # yes / no
        'UnitCells': '1 1 1',                                           # int int int
        'GasComposition': {'CO2': 1.0},                                 # dict
        'SpacingVDWGrid': 0.1,                                          # float
        'SpacingCoulombGrid': 0.1,                                      # float
        'UseTabularGrid': 'no',                                         # yes / no
        'NumberOfGrids': 0,                                             # int
        'GridTypes': '',                                                # string
        'Movies': 'no',                                                 # yes / no
        'WriteMoviesEvery': 0,                                          # int
        'ComputeDensityProfile3DVTKGrid': 'no',                         # yes / no
        'DensityProfile3DVTKGridPoints': '100 100 100',                 # int int int
        'WriteDensityProfile3DVTKGridEvery': 100,                       # int
        'RestartFile': 'no',                                            # yes / no
    }

    # Update the dictionary with the kwargs
    CALC_DICT.update(kwargs)

    if 'UnitCells' in kwargs:
        if isinstance(kwargs['UnitCells'], list):
            kwargs['UnitCells'] = ' '.join(map(str, kwargs['UnitCells']))
        if isinstance(kwargs['UnitCells'], int):
            kwargs['UnitCells'] = ' '.join(map(str, [kwargs['UnitCells']] * 3))
        if isinstance(kwargs['UnitCells'], str) and kwargs['UnitCells'].lower() == 'auto':

            maxCutOff = max([CALC_DICT['CutOffVDW'],
                             CALC_DICT['CutOffChargeCharge'],
                             CALC_DICT['CutOffChargeBondDipole'],
                             CALC_DICT['CutOffBondDipoleBondDipole']])

            CALC_DICT['UnitCells'] = calculate_UnitCells(
                os.path.join(path, FrameworkName.rstrip('.cif') + '.cif'),
                maxCutOff)

    if 'ExternalPressure' in kwargs:
        if isinstance(kwargs['ExternalPressure'], list):
            CALC_DICT['ExternalPressure'] = ' '.join(map(str, kwargs['ExternalPressure']))

        elif isinstance(kwargs['ExternalPressure'], int):
            CALC_DICT['ExternalPressure'] = float(kwargs['ExternalPressure'])

        elif isinstance(kwargs['ExternalPressure'], float):
            CALC_DICT['ExternalPressure'] = kwargs['ExternalPressure']

        elif isinstance(kwargs['ExternalPressure'], str):
            CALC_DICT['ExternalPressure'] = ' '.join(kwargs['ExternalPressure'].split(','))

    # Create file header as string
    GCMC_InputFile = dedent("""\
    SimulationType                      MonteCarlo
    NumberOfCycles                      {NumberOfCycles}
    NumberOfInitializationCycles        {NumberOfInitializationCycles}
    PrintEvery                          {PrintEvery}
    PrintPropertiesEvery                {PrintPropertiesEvery}

    RestartFile                         {RestartFile}

    ForceField                          {ForceField}
    CutOffVDW                           {CutOffVDW}
    CutOffChargeCharge                  {CutOffChargeCharge}
    CutOffChargeBondDipole              {CutOffChargeBondDipole}
    CutOffBondDipoleBondDipole          {CutOffBondDipoleBondDipole}
    ChargeMethod                        Ewald
    EwaldPrecision                      {EwaldPrecision}

    Framework                           0
    FrameworkName                       {FrameworkName}
    HeliumVoidFraction                  {HeliumVoidFraction}
    ExternalTemperature                 {ExternalTemperature}
    ExternalPressure                    {ExternalPressure}
    UseChargesFromCIFFile               {UseChargesFromCIFFile}
    UnitCells                           {UnitCells}

    """).format(**CALC_DICT)

    if CALC_DICT['UseTabularGrid'] == 'yes':

        # Get the pseudoatoms number and types
        for gas in list(CALC_DICT['GasComposition'].keys()):
            pseudo_atoms = get_pseudoatoms(gas)
            CALC_DICT['NumberOfGrids'] += len(pseudo_atoms)           # int
            CALC_DICT['GridTypes'] += ' '.join(pseudo_atoms) + ' '    # string

        GCMC_InputFile += dedent("""\
        NumberOfGrids                       {NumberOfGrids}
        GridTypes                           {GridTypes}
        SpacingVDWGrid                      {SpacingVDWGrid}
        SpacingCoulombGrid                  {SpacingCoulombGrid}
        UseTabularGrid                      {UseTabularGrid}

        """).format(**CALC_DICT)

    if CALC_DICT['ComputeDensityProfile3DVTKGrid'] == 'yes':
        GCMC_InputFile += dedent("""\
        ComputeDensityProfile3DVTKGrid      yes
        DensityProfile3DVTKGridPoints       {DensityProfile3DVTKGridPoints}
        WriteDensityProfile3DVTKGridEvery   {WriteDensityProfile3DVTKGridEvery}

        """).format(**CALC_DICT)

    if CALC_DICT['Movies'] == 'yes':
        GCMC_InputFile += dedent("""\
        Movies yes
        WriteMoviesEvery       {WriteMoviesEvery}

        """).format(**CALC_DICT)

    # Create component list as string
    for name, fraction in CALC_DICT['GasComposition'].items():

        number_of_components = len(CALC_DICT['GasComposition'])
        index_of_component = list(CALC_DICT['GasComposition']).index(name)

        # Append component string block to input file
        GCMC_InputFile += dedent(f"""\
        Component {index_of_component} MoleculeName                  {name}
                    MolFraction                  {fraction}
                    MoleculeDefinition           TraPPE
                    SwapProbability              0.5
                    TranslationProbability       0.3
                    RotationProbability          0.2
                    IdentityChangeProbability    0.1
                        NumberOfIdentityChanges    {number_of_components}
                        IdentityChangesList        {' '.join(map(str, range(number_of_components)))}
                    CreateNumberOfMolecules      0

        """)

    # Write input to file
    with open(os.path.join(path,
                           "simulation.input"), 'w') as f:
        f.write(GCMC_InputFile)
