import numpy as np
from scipy.optimize import curve_fit
from statsmodels.tsa.stattools import adfuller


def exp_decay(t, tau):
    """
    Simple function to model a exponential decay.

    Parameters
    ----------
    t : array
        Time data
    tau : float
        Decay rate of the exponential
    Returns
    -------
    Exponential curve as a NumPy array.
    """

    return np.exp(-t/tau)


def check_consistency(data):
    """
    Checks the consistency of the input data.

    Parameters
    ----------
    data : array
        Array with the data
    Returns
    -------
    consistent : bool
        Boolean indicating if the input data is consistent
    data_array : array
        NumPy array containing the data in the correct format for the next steps.
    """

    # Try to convert the data to a NumPy Array
    try:
        data_array = np.array(data).astype(float)
        # Remove wrong nested lists
        data_array = data_array.squeeze()
    except ValueError:
        print('Input data must be an array of float numbers!')
        print('The following data was passed:')
        print(data)

        # Replace the incorrect data with a array with a NaN value
        data_array = np.array(np.nan)

    # Check if input data is unidimentional
    if data_array.ndim != 1:
        raise Exception(f'Input data must be 1D. {data_array.ndim}D data used instead!')

    # Check if all the data is finite
    is_all_finite = np.all(np.isfinite(data_array))

    # Check if the data is not an array filled with zeros
    is_all_zero = np.all((data_array == 0))

    return is_all_finite, is_all_zero, data_array


def batch_average_data(data, batch_size=1):
    """
    Converts the data to batch averages with a given batch size.

    Parameters
    ----------
    data : array
        Array with the data
    batch_size : int
        Size of the batch to take the averages
    Returns
    -------
    averaged_batches : array
        Array containig the batch-averaged data
    """

    if batch_size > 1:
        # Trucate the data to allow a closed batch.
        # Be aware that this will remove the last points to make a closed batch
        truncated_data = data[:int(
            np.floor(len(data) / batch_size) * batch_size)]

        # Reshape the data to create batch of size m.
        reshaped_data = np.reshape(truncated_data, (-1, batch_size))

        # Get the average of each batch
        averaged_batches = np.array([np.average(i) for i in reshaped_data])

        return averaged_batches

    else:
        return data


def calculate_MSEm(data, batch_size=1):
    """
    Calculates the m-Marginal Standard Error (MSEm) for a simulation data
    with batch size equals to m. m=1 reduces to the original MSER.

    Parameters
    ----------
    data : array
        Array with the data
    batch_size : int
        Size of the batch to take the averages
    Returns
    -------
    MSE : array
        Array containig the Marginal Standard Error data
    """

    # Convert data to n-blocked average
    block_average = batch_average_data(data, batch_size)

    # Get the size of the data
    n = len(block_average)

    # Creates a empty list to store the MSE values
    MSE = []

    # Iterate over data index and calculates the average from k to n-2
    for k in range(n - 2):
        # Truncate data on k and convert it to a numpy array
        truncated_data = np.array(block_average[k:])

        # Get the average of the truncated data
        Y_nk = np.average(truncated_data)

        # Calculates the sum of the squared diference
        sum_sq_diff = np.sum([(j - Y_nk)**2 for j in truncated_data])

        # Calculate the k-th Marginal standard error
        g_k = sum_sq_diff / (n - k)**2

        # Add the k-th to MSE array
        MSE += [g_k]

    return np.array(MSE)


def MSERm_index(MSEm, batch_size=1):
    """
    Applies the m-Marginal Standard Error Rule (MSERm) to the SERm data to get
    the position where equilibrated data starts.

    Parameters
    ----------
    MSEm : array
        Marginal Standard Error applied to the data
    batch_size : int
        Size of the batch to take the average
    Returns
    -------
    equilibrated_index : int
        Index of the start of equilibrated data
    """
    # Remove potential too low values that apears artificially on last points
    MSEm = np.where(np.array(MSEm) < 1e-9,   # where value < 1e-9
                    max(MSEm),               # replace for max(MSEm)
                    np.array(MSEm))          # on MSEm array

    equilibrated_index = np.argmin(MSEm)*batch_size

    return equilibrated_index


def MSERm_LLM_index(MSEm, batch_size=1):
    """
    Applies the LLM version of m-Marginal Standard Error Rule (MSERm) to the SERm
    data to get the position where equilibrated data starts. This method gets
    the first minimum on Marginal Standard Error curve and assumes it is the
    start of equilibriation. It is a better option for complicated adsorptions
    like water close to condensation.

    Parameters
    ----------
    MSEm : array
        Marginal Standard Error applied to the data
    batch_size : int
        Size of the batch to take the average
    Returns
    -------
    t0 : int
        Start of the LLM equilibrated data
    """

    # Search for the first mininum on the MSEm data
    i = 0
    while MSEm[i+1] < MSEm[i]:
        i += 1
    # Correct for the batch size
    t0 = i*batch_size

    return t0


def calc_equilibrated_average(data, eq_index, uncertainty='uSD', ac_time=1):
    '''
    Calculates the average and uncertainty on the equilibrated part
    of the data.

    Parameters
    ----------
    data : array
        Array with the data
    eq_index : int
        Index of the start of equilibrated data.
    uncertainty : str
        String for selecting Standard Error (SE), Standard Deviation (SD), or its
        uncorrelated versions uSD and uSE as the default uncertainty of the average.
    ac_time : int
        Autocorrelation time
    Returns
    -------
    equilibrated_average : float
        Average on the equilibrated data
    equilibrated_uncertainty : float
        Uncertainty of the average calculation
    '''

    if uncertainty not in ['SD', 'SE', 'uSD', 'uSE']:
        raise Exception(f"""{uncertainty} is not a valid option!
            Only Standard Deviation (SD), Standard Error (SE), uncorrelated
            Standard Deviation (uSD), and uncorrelated Standard Error (uSE)
            are valid options.""")

    # Remove the initial transient of the data
    equilibrated_data = data[eq_index:]

    # Calculates the average on the equilibrated data
    equilibrated_average = np.average(equilibrated_data)

    # Calculate the standad deviation on the equilibrated data
    if uncertainty == 'SD':
        equilibrated_uncertainty = np.std(equilibrated_data)

    # Calculate the Standard Error
    elif uncertainty == 'SE':
        equilibrated_uncertainty = np.std(equilibrated_data) / np.sqrt(len(equilibrated_data))

    # Calculate the uncorrelated Standard Error
    elif uncertainty == 'uSD':
        # Divide the equilibrated_data on uncorrelated chunks
        uncorr_batches = batch_average_data(equilibrated_data,
                                            np.ceil(ac_time).astype(int))

        # Calculate the standard deviation on the uncorrelated chunks
        equilibrated_uncertainty = np.std(uncorr_batches)

    # Calculate the uncorrelated Standard Error
    elif uncertainty == 'uSE':
        # Divide the equilibrated_data on uncorrelated chunks
        uncorr_batches = batch_average_data(equilibrated_data,
                                            np.ceil(ac_time).astype(int))

        # Calculate the standard error of the mean on the uncorrelated chunks
        equilibrated_uncertainty = np.std(uncorr_batches) / np.sqrt(len(uncorr_batches))

    return equilibrated_average, equilibrated_uncertainty


def calc_autocorrelation_time(data):
    """
    Calculates the autocorrelation time of a equilibrated data.
    Autocorrelation is expected to fall off exponentially at long times

    Parameters
    ----------
    data : array
        Array of data to calculate the integrated autocorrelation time
    Returns
    -------
    autocorrelation_time : float
        Autocorrelation time
    uncorrelated_samples : float
        Number of uncorrelated samples
    """

    # Check the consistency of the time_serie
    is_all_finite, is_all_zero, data_array = check_consistency(data)

    if is_all_finite is False or is_all_zero is True:
        return 0, 0

    try:
        # Calculates the ACF using numpy
        data_std = data_array - np.mean(data_array)
        data_norm = np.sum(data_std ** 2)
        ACF = np.correlate(data_std, data_std, mode='full')/data_norm
        ACF = ACF[int(ACF.size/2):]

        # Fit a exponential decay to ACF
        x = np.arange(len(ACF))
        [tau], _ = curve_fit(exp_decay,  x,  ACF)

        # Calculate autocorrelation time as the half-live of ACF exponential decay
        autocorrelation_time = np.ceil(tau*np.log(2))

    except (RuntimeError, ValueError) as Error:
        # If the if the least-squares minimization fails, set the autocorrelation_time to 1.
        # This can happen if the ACF data do not present a exponential decay
        autocorrelation_time = 1
        print('The least-squares minimization failed! Please check the data.')
        print(Error)

    # Calculate the number of uncorrelated data
    uncorrelated_samples = data_array.size / autocorrelation_time

    return autocorrelation_time, uncorrelated_samples


def apply_ADF_test(equilibrated_data, verbosity=True):
    """
    Applies the Augmented Dickey-Fuller Test on the equilibrated data

    Parameters
    ----------
    equilibrated_data : array
        Array with the equilibrated data
    verbosity : bool
        Boolean to control the output printing
    Returns
    -------
    ADFTestResults : dict
        Dictionary containg the ADF test results
    output : str
        String containg the output
    """
    adf, p, usedlag, n_obs, cv, icbest = adfuller(equilibrated_data, autolag='AIC')

    ADFTestResults = {'adf': adf,
                      'pvalue': p,
                      'usedlag': usedlag,
                      'n_obs': n_obs,
                      'critical_values': cv,
                      'icbest': icbest}

    output = f"""
                           Augmented Dickey-Fuller Test
==============================================================================
Test statistic for observable: {adf}
P-value for observable: {p}
The number of lags used: {usedlag}
The number of observations used for the ADF regression: {n_obs}
Cutoff Metrics :
"""
    for k, v in cv.items():
        conf = 100 - int(k.rstrip('%'))
        if v < adf:
            output += f"{k:>4}: {v:9.6f} | The data is not stationary with {conf} % confidence\n"
        else:
            output += f"{k:>4}: {v:9.6f} | The data is stationary with {conf} % confidence\n"

    if verbosity:
        print(output)

    return ADFTestResults, output


def equilibrate(input_data,
                LLM=False,
                batch_size=1,
                ADF_test=True,
                uncertainty='uSD',
                print_results=True):
    """
    Wrap function to apply MSER to an input_data array.

    Parameters
    ----------
    input_data : array
        Array with the original data
    LLM : bool
        Boolean to control usage of the LLM variation of MSER
    batch_size : int
        Size of the batch to take the average
    ADF_test : bool
        Boolean to control usage ADF test
    uncertainty : str
        String for selecting Standard Error (SE), Standard Deviation (SD), or its
        uncorrelated versions uSD and uSE as the default uncertainty of the average.
    print_results : bool
        Boolean to control printing of the results
    Returns
    -------
    results_dict : dict
        Dictionary containg resunts of MSER
    """

    # Check the consistency of the time_serie
    is_all_finite, is_all_zero, array_data = check_consistency(input_data)

    # Returns NaN if any of the time_series data is not a finite number
    if is_all_finite is False:
        results_dict = {'MSE': np.nan,
                        't0': np.nan,
                        'average': np.nan,
                        'uncertainty': np.nan,
                        'equilibrated': np.nan,
                        'ac_time': np.nan,
                        'uncorr_samples': np.nan}
        return results_dict

    # Returns zero if all the data in time_series is zero
    if is_all_zero:
        results_dict = {'MSE': np.zeros(len(array_data)),
                        't0': 0,
                        'average': 0,
                        'uncertainty': 0,
                        'equilibrated': np.zeros(len(array_data)),
                        'ac_time': 0,
                        'uncorr_samples': 0}
        return results_dict

    # Check if the input parameters are what is expected
    assert isinstance(LLM, bool), 'LLM should be True or False'
    assert isinstance(batch_size, int), 'batch_size should be an int'
    assert isinstance(ADF_test, bool), 'ADF_test should be True or False'
    assert isinstance(print_results, bool), 'print_results should be True or False'

    # Check if the uncertainty is a valid option
    if uncertainty not in ['SD', 'SE', 'uSD', 'uSE']:
        raise Exception(f"""{uncertainty} is not a valid option!
            Only Standard Deviation (SD), Standard Error (SE), uncorrelated
            Standard Deviation (uSD), and uncorrelated Standard Error (uSE)
            are valid options.""")

    # Calculate the Marginal Standard Error curve
    MSEm_curve = calculate_MSEm(array_data, batch_size=batch_size)

    if LLM is False:
        # Apply the MSER to get the index of the start of equilibrated data
        t0 = MSERm_index(MSEm_curve, batch_size=batch_size)

    if LLM is True:
        # Apply the MSER-LLM to get the index of the start of equilibrated data
        t0 = MSERm_LLM_index(MSEm_curve, batch_size=batch_size)

    # Calculate autocorrelation time and the number of uncorrelated samples
    equilibrated = array_data[t0:]
    ac_time, uncorr_samples = calc_autocorrelation_time(equilibrated)

    # Calculates the average and standard deviation on the equilibrated data
    average, avg_uncertainty = calc_equilibrated_average(array_data,
                                                         t0,
                                                         uncertainty,
                                                         ac_time)

    # Create a dictionary with the results
    results_dict = {'MSE': MSEm_curve,
                    't0': t0,
                    'average': average,
                    'uncertainty': avg_uncertainty,
                    'equilibrated': equilibrated,
                    'ac_time': ac_time,
                    'uncorr_samples': uncorr_samples}

    eq_ratio = 100 * (len(array_data) - t0) / len(array_data)

    if print_results:
        print(f"""                            pyMSER Equilibration Results
==============================================================================
Start of equilibrated data:          {t0} of {len(array_data)}
Total equilibrated steps:            {len(array_data) - t0}  ({eq_ratio:.2f}%)
Average over equilibrated data:      {average:.4f} ± {avg_uncertainty:.4f}
Number of uncorrelated samples:      {uncorr_samples:.1f}
Autocorrelation time:                {ac_time:.1f}
==============================================================================""")

    if ADF_test:
        # Apply the Augmented Dickey-Fuller Test on the equilibrated data
        ADFTestResults, output_text = apply_ADF_test(equilibrated,
                                                     verbosity=print_results)
        results_dict.update(ADFTestResults)

    return results_dict
