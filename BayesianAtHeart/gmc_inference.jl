################################################################################
# Bayesian estimation of heart rate dynamics
#
# This script runs the estimation of heart rate trajectories.  It is meant to
# be used in conjunction with the script `generate_hr.py', from which it is
# called.
#
# The script relies primarily in the open source package PointProcessInference,
# which was introduced in the following article:
#
# Gugushvili, S., van der Meulen, F., Schauer, M., & Spreij, P. (2018).  Fast
# and scalable non-parametric Bayesian inference for Poisson point processes.
# arXiv preprint arXiv:1804.03616.
#
# Fernando Rosas & Pedro Mediano, March 2023
#
# ##############################################################################

using PointProcessInference 
const PPI = PointProcessInference
using Random
using Statistics
using CSV
using DataFrames
using Distributions

function generate_hr(data, IT, tau, theta, Nr, Nd, verbose)
  # Building a prior from basic data stats
  RR = diff(data);
  ss = fit_mle(Gamma,1 ./ RR);
  α1 = αind = ss.α;
  β1 = βind = 1/ss.θ;

  # Setting the sampling frequency (Fs=3Hz)
  N = round(Int64, (maximum(data) - minimum(data))*3 );

  # Create a conteiner for the IT sampled trajectories, each of length N
  M = zeros(N,IT+1);

  # Sample each of the IT trajectories
  for it=1:IT

    if verbose>0
      print("Iteration n" * string(it) * "\n");
    end

    # Main function call to GMC inference routine
    # data: interbeat sequence to use as data for the algorithm
    # T0: starting time
    # N: number of temporal bins
    # title: name of the output column
    # samples: indices of the realisations of the Gibbs sampler to place in the output
    # Nr: number of runs of Gibbs sampler
    # α1, β1: hyperparameters of prior of timepoint t0
    # αind, βind: hyperparameters of the prior of the whole Markov chain
    # τ: value of parameter of Metropolis-Hastings procedure, related to the prior of gamma
    # Π: prior of Gamma
    F = PPI.inference(data; T0=minimum(data), N=N, title = "HR", samples=1:1:Nr, α1=α1, β1=β1, αind=αind, βind=βind, τ=tau, Π=Exponential(theta), verbose=false);

    X = F.ψ; # extract the obtained values 
    X = X[Nd:end,:]; # discard the first Nd realisations
    M[:,it+1] = mean(X, dims=1) * 60; # take the mean value of the remaining realisations

    if it==1 # make a copy of the resulting time indices (which are the same for all iterations)
      t_index = collect(F.breaks);
      t_mean = (t_index[1:end-1]+t_index[2:end])/2
      M[:,it] = t_mean;
    end

  end

  return M

end 


function generate_hr_from_files(files, IT, tau, theta, Nr, Nd, verbose)
  """
  Function to generate HR trajectories using a gamma Markov chain (GMC) model

  Parameters
  ----------
  files: list
    List of strings specifying the files of data that will be used as inputs
  IT  : int
          Number of sampled trajectories extracted via the Gibbs sampler
      tau : int
          Random walk step-size in Metropolis-Hastings step for estimating gamma
  theta: float
          Hyperparameter of prior of gamma. In general, theta=0.01 gives high bandwidth, theta=10 gives low bandwidth
      Nr  : int
          Number of runs of the Gibbs sampler per estimated trajectory
  Nd  : int
          Number of runs of the Gibbs sampler discarded before calculating average

  Returns
  -------
  output : list
    List of arrays of IT+1 columns, where the first is the time index and the other IT are the estimated trajectories
  """
  Random.seed!(1234); # set a fixed seed for reproducibility
  output = [];
  for f in files

    if verbose>0
      full_name = split(f,'/')[end];
      name = split(full_name,'.')[1];
      print("\nRunning Bayesian estimation of " * name * "\n");
    end

    # Loading the corresponding data
    d = CSV.read( f, DataFrame; header=false);
    data = d[:,1]
    M = generate_hr(data, IT, tau, theta, Nr, Nd, verbose)
    # Append M to output list
    push!(output, M)

    #saving the generated estimation to csv file
    df = DataFrame(M)
    new_name = string(split(f,".")[1],"_hr_estimation.csv")
    CSV.write(new_name,df,header=false)
    println("Estimation saved in ", new_name)

  end 

  return output

end


# Both executive functions read .CSV files containg R_peaks data and save
# .CSV files containg HR trajectories using <original path>_hr_estimation.csv pattern

function run()
    files = [
        # "data/R_peaks/VRCC_S01.csv",
        # "data/R_peaks/VRCC_S02.csv",
        # "data/R_peaks/VRCC_S03.csv",
        # "data/R_peaks/VRCC_S04.csv",
        # "data/R_peaks/VRCC_S05.csv",
        # "data/R_peaks/VRCC_S06.csv",
        # "data/R_peaks/VRCC_S07.csv",
        # "data/R_peaks/VRCC_S08.csv",
        # "data/R_peaks/VRCC_S09.csv",
        # "data/R_peaks/VRCC_S10.csv",
  
        "data/R_peaks/1_ECG_Misia_baseline_r_peaks.csv",
        "data/R_peaks/1_ECG_Misia_main_r_peaks.csv",
        "data/R_peaks/2_ECG_Michal_baseline_r_peaks.csv",
        "data/R_peaks/2_ECG_Michal_main_r_peaks.csv"    ]

    outs = generate_hr_from_files(files,3,1,1,20000,5000,1)
end

function run_example()
    outs = generate_hr_from_files(["data/example/data.csv"],3,1,1,20000,5000,1)
end
