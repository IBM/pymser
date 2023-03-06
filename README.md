# pyMSER

A Python library to apply the [Marginal Standard Error Rule (MSER)](https://doi.org/10.1177/003754979706900601) for transient regime detection and truncation on Grand Canonical Monte Carlo adsorption simulations.

## Dependencies

* [NumPy](https://numpy.org) is the fundamental package for scientific computing with Python.
* [SciPy](https://scipy.org/) is a collection of fundamental algorithms for scientific computing in Python.
* [statsmodels](https://www.statsmodels.org/) is a python module that provides classes and functions for the estimation of many different statistical models, as well as for conducting statistical tests, and statistical data exploration.

## Developer tips

These tips are not mandatory, but they are a sure way of helping you develop the code while maintaining consistency with the current style, structure and formatting choices.

### Coding style guide

We recommend these tools to ensure code style compatibility.

* [autopep8](https://pypi.org/project/autopep8/) automatically formats Python code to conform to the PEP8 style guide.
* [Flake8](https://flake8.pycqa.org) is your tool for style guide enforcement.

## Installation

### Option 1: Using `setup.py`

Clone `pymser` repository if you haven't done it yet.

Go to `pymser`'s root folder, there you will find `setup.py` file, and run the command below:

```Shell
python setup.py install
```

### Option 2: Using pip/pipenv to install from Pypi.org

If you intend to use `pipenv`, please add the following to your `Pipfile`:

```Pipfile
[[source]]
url = "https://pypi.org/simple"
verify_ssl = true
name = "pypi"

[packages]
pymser = "*"
```

If you intend to use `pip`, please run the command below:

```Shell
pip install pymser
```

### Option 3: Using pip to install directly from the GitHub repo

You can run

```Shell
pip install git+https://github.com/IBM/pymser.git
```

and then you will be prompted to enter your GitHub username and password/access token.

If you already have a SSH key configured, you can run

```Shell
pip install git+ssh://git@github.com/IBM/pymser.git
```

### Option 4: Using pip/pipenv to install from Artifactory

Log into Artifactory and access your user profile. There you will find your API key and username. Then export your credentials as environment variables for later use in the installation process.

```Shell
export ARTIFACTORY_USERNAME=username@email.com
export ARTIFACTORY_ACCESS_TOKEN=your-access-token
export ARTIFACTORY_URL=your-artifactory-url
```

If you intend to use `pipenv`, please add the following to your `Pipfile`:

```Pipfile
[[source]]
url = "https://$ARTIFACTORY_USERNAME:$ARTIFACTORY_ACCESS_TOKEN@$ARTIFACTORY_URL"
verify_ssl = true
name = "artifactory"

[packages]
pymser = {version="*", index="artifactory"}
```

If you intend to use `pip`, please run the command below:

```Shell
pip install pymser --extra-index-url=https://$ARTIFACTORY_USERNAME:$ARTIFACTORY_ACCESS_TOKEN@$ARTIFACTORY_URL
```

## Usage example

This is a small example of how to use our package:

```Python
>>> import pymser
>>> import pandas as pd
>>>
>>> # Reads the example file using pandas
>>> df = pd.read_csv('example_data/Cu-BTT_500165.0_198.000000.csv')
>>>
>>> # Apply the MSER to get the index of the start of equilibrated data
>>> results = pymser.equilibrate(df['mol/kg'], LLM=False, batch_size=1, ADF_test=True, uncertainty='uSD', print_results=True)

                            pyMSER Equilibration Results
==============================================================================
Start of equilibrated data:          13368 of 48613
Total equilibrated steps:            35245  (72.50%)
Equilibrated:                        Yes
Average over equilibrated data:      22.4197 ± 0.1905
Number of uncorrelated samples:      22.3
Autocorrelation time:                1579.0
==============================================================================

                           Augmented Dickey-Fuller Test
==============================================================================
Test statistic for observable: -3.926148246630434
P-value for observable: 0.001850619485090052
The number of lags used: 46
The number of observations used for the ADF regression: 35198
Cutoff Metrics :
  1%: -3.430536 | The data is stationary with 99 % confidence
  5%: -2.861622 | The data is stationary with 95 % confidence
 10%: -2.566814 | The data is stationary with 90 % confidence
```

You can also access our [tutorial](pymser_tutorial.ipynb).

## Python package deployment

### Deploying to Artifactory

We have an automated CI/CD pipeline running on TravisCI that takes every single `git push` event and executes the build/test/deploy instructions in the `.travis.yml`. If you are deploying `master` or `release` branches, a Python package will be generated and published to a private Pypi registry on Artifactory.

### Deploying to Pypi

We have an automated CI/CD pipeline running on TravisCI that takes every single `git push` event and executes the build/test/deploy instructions in the `.travis.yml`. If you are deploying `main` branch, a Python package will be generated and published to Pypi.org registry.
