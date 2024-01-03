# Changelog

All notable changes to this project will be documented in this file.

## v1.0.20

- Use `pyTorch` to calculate the MSE curve reducing the overall time of the calculation by approximately one order of magnitude on large arrays
- Use scipy to calculate the autocorrelation time reducing the overall time of the calculation by approximately two order of magnitudes on large arrays

## v1.0.18

- Use `nanmean` and `nanstd` instead of `mean` and `std` to avoid erros when there are `NaN` values in the data
- Add equilibration status on the print report and a warning if the equilibration is not reached

## v1.0.8

- Downgrade requirements for Python from python>=3.10 to python>=3.9
- Add the Standard Error (SE) as possible uncertainty of the average
- Add the uncorrelated Standard Error (uSE) as possible uncertainty of the average
- Add the uncorrelated Standard Deviation (uSD) as possible uncertainty of the average and set it as default
- Small bug fixes

## v1.0.2

- Add files to github repository
- Prepare for Pypi.org release
