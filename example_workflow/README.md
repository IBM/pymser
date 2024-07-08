# Tutorial: Using pyMSER with RASPA2 to Run Fixed Number of Production Cycles

## Introduction

Welcome to this tutorial on using pyMSER with RASPA2 to run a fixed number of production cycles after automatically detecting the equilibrated portion of the simulation. This guide assumes that you have a basic understanding of Python and Grand Canonical Monte Carlo (GCMC) simulations using RASPA. The aim is to streamline simulation workflows by leveraging the capabilities of pyMSER for equilibration detection, ensuring accurate and efficient production runs.

## Installation

To get started, ensure that your environment is set up with the necessary dependencies. The provided `environment.yml` file includes all the required packages:

* [NumPy](https://numpy.org) is the fundamental package for scientific computing with Python.
* [Pandas](https://pandas.pydata.org) is a fast, powerful, flexible, and easy-to-use open-source data analysis and data manipulation library built on top of NumPy.
* [Gemmi](https://gemmi.readthedocs.io/en/latest/) is a Python library for handling macromolecular structures.
* [pyMSER](https://pypi.org/project/pymser/) is a Python library to apply the Marginal Standard Error Rule (MSER) for transient regime detection and truncation on GCMC adsorption simulations.
* [RASPA2](https://pypi.org/project/RASPA2/) is a Python interface to the RASPA2 molecular simulation package.

To create the environment, run the following command:

```sh
conda env create -f environment.yml
```

Activate the environment with:

```sh
conda activate pymser
```

## How to Run

To run a fixed number of production cycles using pyMSER with RASPA2, you can use the provided `run.py` script. This script automates the process of equilibration detection and production run based on the MSER rule. The script takes the framework name, external pressure, and simulation type as input arguments.

```sh
python run.py --FrameworkName 'MgMOF-74' --ExternalPressure 1e4 --NumberOfProdCycles 1500 --AddCycles 1000 --GasComposition '{"CO2":0.5,"N2":0.5}' 'GCMC''
```

The script will perform the following steps:

1. Create the necessary directories and input files for the simulation.
2. Run RASPA for `AddCycles` cycles.
3. Apply the MSER rule to detect the equilibrated portion of the simulation.
4. If the desired number of production cycles is not reached, run additional `AddCycles cycles until the target is achieved.
5. Parse the output files and save the results.

The output should be:

```sh
==============================================================================
                    Automatic GCMC Simulation with pyMSER
==============================================================================
Framework Name        : MgMOF-74
External Temperature  : 298.0 K
External Pressure     : 10000.0 Pa
Gas Composition       : {'CO2': 0.5, 'N2': 0.5}
Desired Prod. Cycles  : 1500
Output Path           : GCMC
==============================================================================

    Running RASPA simulation...
    
       > Running iteration 1...
       > Found only 999/1500 production cycles. Running more 1000 cycles.
       > Running iteration 2...
       > Success! Found 1999 production cycles. Analyzing final data...


                            pyMSER Equilibration Results
==============================================================================
Start of equilibrated data:          1 of 2000
Total equilibrated steps:            1999  (99.97%)
Equilibrated:                        Yes
Average over equilibrated data:      3.5916 ± 1.9509
Number of uncorrelated samples:      1999.0
Autocorrelation time:                1.0
==============================================================================

                           Augmented Dickey-Fuller Test
==============================================================================
Test statistic for observable: -40.707108860029564
P-value for observable: 0.0
The number of lags used: 0
The number of observations used for the ADF regression: 3998
Cutoff Metrics :
  1%: -3.431987 | The data is stationary with 99 % confidence
  5%: -2.862263 | The data is stationary with 95 % confidence
 10%: -2.567155 | The data is stationary with 90 % confidence

==============================================================================
Component 0 [CO2]
---------------------------------------------------------------------------
Average loading absolute [molecules/unit cell]         0.1947205551 +/-          0.1135451716
Average loading absolute [mol/kg framework]            0.0891424690 +/-          0.0519806290
Average loading absolute [mg/g framework]              3.9231154891 +/-          2.2876414927
Average loading absolute [cm^3 STP/gr]                 1.9980371381 +/-          1.1650925581
Average loading absolute [cm^3 STP/cm^3]               1.7647568024 +/-          1.0290624625
Enthalpy of adsorption [kJ/mol]                      -17.8534776591 +/-          0.3852622588
==============================================================================
Component 1 [N2]
---------------------------------------------------------------------------
Average loading absolute [molecules/unit cell]         0.0297574394 +/-          0.0428471501
Average loading absolute [mol/kg framework]            0.0136228639 +/-          0.0196152930
Average loading absolute [mg/g framework]              0.3816227347 +/-          0.5494910502
Average loading absolute [cm^3 STP/gr]                 0.3053425404 +/-          0.4396567026
Average loading absolute [cm^3 STP/cm^3]               0.2696923470 +/-          0.3883246922
Enthalpy of adsorption [kJ/mol]                      -11.9135661990 +/-          3.2649508929
==============================================================================
```

### Command Line Options

You have several options to control the simulation parameters. Below is a detailed explanation of each option available in the `run.py` script:

- **output_folder (required)**
  - Type: `str`
  - Help: Directory to save the files of the calculations. This directory should contain the `cif` file of the framework and the force field files.

- **--FrameworkName (required)**
  - Type: `str`
  - Help: Name of the framework to be simulated

- **--ExternalPressure (required)**
  - Type: `float`
  - Help: External pressure in Pascal

#### Optional Parameters

- **--NumberOfProdCycles**
  - Type: `int`
  - Default: `5000`
  - Help: Number of desired production cycles

- **--AddCycles**
  - Type: `int`
  - Default: `1000`
  - Help: Number of additional tentative cycles if the desired number of production cycles has not been achieved

- **--ExternalTemperature**
  - Type: `float`
  - Default: `298.0`
  - Help: Temperature of the simulation in Kelvin

- **--UnitCells**
  - Type: `str`
  - Default: `auto`
  - Help: Number of unit cells to be simulated. Can be "auto" or a string of comma-separated values. E.g., "3,3,1"

- **--GasComposition**
  - Type: `str`
  - Default: `{"CO2": 1.0}`
  - Help: Type of gas composition for the simulation as a dictionary. E.g., '{"CO2": 0.5, "N2": 0.5}'

- **--UseChargesFromCIFFile**
  - Help: Use charges from CIF file.


With these options, you can customize your simulation to fit your specific needs, ensuring an efficient and accurate GCMC simulation workflow.
