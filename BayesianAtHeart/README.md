# HR ENTROPY ESTIMATION 
This code allows to estimate HR dynamics following a Bayesian paradigm 
## Download, installation, and code examples

Steps to calculate HR_Entropy
1) Convert ACQ file(s) to CSV and calculate cumsum to get R_peaks using `R_peaks_pure.ipynb`
2) Run `gmc_inference.jl` which runs on Julia. 

```julia (in Repl)
    using Pkg
    Pkg.activate('.')
    Pkg.instantiate()
    include("gmc_inference.jl")
    run()
```
3) Calculate Bayesian (and Frequentiest) estimation of heart rate dynamics using `generate_hr.py` file
4) Calculate HR Entropy with `calculate_HR_entropy.py` file


Note: you might have to play around with Java a bit to execute point 4)