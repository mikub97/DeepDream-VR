"""
Bayesian estimation of heart rate dynamics

This script takes files with sequences of inter-beat intervals as inputs, and
generate heart rate trajectories that are likely to have generated such data.

Fernando Rosas & Pedro Mediano, March 2023
"""

#######################################################
# Loading packages
#######################################################
import pandas as pd
import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import seaborn as sns

#######################################################
# Functions
#######################################################
def hr_interpolate(L_data, fs=1):
    '''
    Function for align data via interpolation

    Parameters
    ----------
    L_data : list of pd.DataFrame
        List of dataframes where each column is a HR timeseries, and index is the time axis
    fs : float
        Sampling frequency of output after interpolation

    Returns
    -------
    L_out : list of pd.DataFrame
        List of dataframes containing HR estimation under equally sampled data        
    '''
    L_out = []
    for k in range(len(L_data)):
        data = L_data[k]
        T = data.index

        # Create common time index
        time_ix = np.arange( np.ceil(T[0]), np.floor(T[-1]), 1/fs)
        df_out = pd.DataFrame( index=time_ix, columns=data.columns)
        df_out.index.name='Time'

        # interpolate data from each sample
        for c in data.columns:
            D = data[c].values
            F = interpolate.interp1d(x=T, y=D, kind='cubic') # preparation
            df_out.loc[time_ix,c] = F(time_ix) # carrying out the actual interpolation
        L_out.append(df_out)

    return L_out


def frequentist_hr(filenames, verbose=1):
    '''
    Function to generate HR trajectory via the standard frequentist approach, 
    using the formula hr = 60 / inter-beat interal

    Parameters
    ----------
    filenames: list of str
        List of strings specifying the files of data that will be used as inputs

    Returns
    -------
    hr_aligned : list of pd.DataFrame
        List of dataframes containing the frequentist HR estimation under equally sampled data
    '''
    L_data = []

    for f in filenames:

        if verbose>0:
            name = f.split('/')[-1].split('.')[0]
            print('\nRunning frequentist estimation of ',name)

        # Load and transform
        S = pd.read_csv(f,header=None, sep=',')[0] # Read the data
        rr = S.diff().iloc[1:] # Calculate inter-beat intervals, discard the first because is NaN
        if rr.min()<0:
            raise ValueError('No useful data found in file '+f)

        HR = 60/rr.values # Calculate the instantaneous HR
        time = S.values[:-1] + rr.values/2 # Attribute the HR value to the midpoint between the beats 
        df = pd.DataFrame(data=HR, index=time)
        df.index.name='Time'
        L_data.append(df)

    # Interpolate data to get regular timepoints with sampling frequency fs
    hr_aligned = hr_interpolate(L_data)

    return hr_aligned


def bayesian_hr(filenames, IT=3, theta=1., tau=1., Nr=20000, Nd=5000, rol=9, w_type='triang', dec=3, fs=1, verbose=1, script_location=None):
    """
    Function to generate Bayesian estimates of HR.

    Parameters
    ----------
    filenames: list of str
        List of strings specifying the files of data that will be used as inputs
    IT  : int
        Number of sampled trajectories extracted via the Gibbs sampler
    theta: float
        Hyperparameter of prior of gamma. In general, theta=0.01 gives high bandwidth, theta=10 gives low bandwidth
    tau : int
        Random walk step-size in Metropolis-Hastings step for estimating gamma
    Nr  : int
        Number of runs of the Gibbs sampler per estimated trajectory
    Nd  : int
        Number of runs of the Gibbs sampler discarded before calculating average
    rol : int
        Number of samples involved in rolling mean used for post-processing smoothing
    w_type: string
        Type of window used for rolling mean
    dec : int
        Downsampling level. Recommended dec=3 for attaining overall sampling of ~1Hz
    fs  : int
        Sampling frequency of final data, after interpolation
    script_location : string
        Location of the file 'gmc_inference.jl'. If none, is assumed that is located in the running folder

    Returns
    -------
    hr_aligned : list of pd.DataFrame
        List of dataframes containing each sampled HR trajectory under equally sampled data
    """
    if script_location is None:
        # If no script_location is provided, the script assumes file is in the running folder
        script_location = 'gmc_inference.jl'

    # Call Julia script
    # generate_hr = J.include(script_location) # Makes a wraper of the Julia function
    # data = generate_hr(filenames, IT, tau, theta, Nr, Nd, verbose) # Calls generate_hr script
    L_data = []

    for f in filenames:
        data = pd.read_csv(f, sep=',', header=None).values

        df = pd.DataFrame(data=data[:,1:], index=data[:,0])
        df.index.name = 'Time'
        df = df.rolling(rol, center=True, win_type=w_type).mean().dropna()
        df = df.iloc[::dec]
        L_data.append(df)
    
    # Interpolate data to get regular timepoints with sampling frequency fs
    hr_aligned= hr_interpolate(L_data, fs)

    return hr_aligned 


###############################
# Example of HR estimation
###############################
if __name__ == '__main__':

    # Define list of filenames

    # patients = [f'VRCC_S0{str(i)}' for i in range(1,10)]+['VRCC_S10']
    # filenames_est = [f'./R_peaks/{p}_hr_estimation.csv' for p in patients]
    # filenames_clean = [f'./R_peaks/{p}.csv' for p in patients]

    filenames = ['1_ECG_Misia_main','2_ECG_Michal_main']

    filenames_est = [f"data/hr_estimations/{f}_r_peaks_hr_estimation.csv" for f in filenames]
    filenames_clean = [f"data/R_peaks/{f}_r_peaks.csv" for f in filenames]

    # Authors example
    # filenames_est = ['./example/data_estimation.csv']
    # filenames_clean = ['./example/data.csv']


    # Calculate Bayesian and frequentist estimation of HR
    bayes_hrs = bayesian_hr(filenames_est, IT=3)
    freq_hrs  = frequentist_hr(filenames_clean)

    assert len(filenames_clean) == len(filenames_est)

    for i in range(len(filenames_est)):
        bayes_hr = bayes_hrs[i]
        freq_hr = freq_hrs[i]

        # Saving results
        print('\nSaving results')
        bayes_hr.to_csv(f'{filenames_est[i].split(".csv")[0]}_bayes_hr.csv')
        freq_hr.to_csv(f'{filenames_est[i].split(".csv")[0]}_freq_hr.csv')

        # Plotting the resulting HR trajectories
        print('\nPloting results')

        sns.set_style('whitegrid')
        sns.set_palette('Set1')

        plt.plot(freq_hr.index, freq_hr.values, linewidth = 1.5, color='black')
        X = bayes_hr.reset_index().melt(id_vars='Time', var_name='Run', value_name='HR')
        g = sns.lineplot(data=X, x='Time', y='HR', errorbar='sd', linewidth=1.5)
        plt.title(filenames_est[i])
        plt.legend(['Frequentist','Bayesian'])
        # plt.show()
        plt.savefig(f"{filenames_est[i].split('.csv')[0]}.png")
        plt.close()

