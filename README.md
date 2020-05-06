# Estimation of BVDV transmission rate using ABC-SMC

### Brief introduction
The purpose of the script is to estimate the within-herd transmission rate (beta) of bovine viral diarrhoea virus (BVDV)
in 9 New Zealand beef herds where animals are extensively grazed. 
We used an approximate Bayesian computation-serial Monte Carlo (ABC-SMC) method to infer the model parameter values.
Detailed description about the ABC-SMC algorhithm can be found in [Toni et al., 2009](https://royalsocietypublishing.org/doi/full/10.1098/rsif.2008.0172).

This script models the spread of BVDV in 9 replacement herds simultaneously. For each herd, x number of persistently infected (PI) animals
is introduced at random time point (tau), and the BVDV Ab seropositivity of 15 random animals is measured 
(1) after the first calving and (2) at pregnancy scanning. The seroconversion status of the 15 animals is compared with the observed
BVDV Ab status from 9 New Zealand beef farms.


### How to run
The simulation model itself was coded in the C programming language, however, the ABC-SMC algorithm was designed to run in R.
In order to run the simulation model in R, the C code should be compiled for [dynamic loading](https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/SHLIB).

For Windows users, open the command prompt at the folder where the c code is located and type;
```sh
R cmd shlib FILENAME.c
```
where "FILENAME.c" is the file name of the C code. To run the command above, one may need to install [R tools](https://cran.r-project.org/bin/windows/Rtools/).


### Notes
The R code was designed to use parallel computation with 10 cores to acheive the computational efficiency. Although one may change the
R code to run the model under a single core computational setting, the running time of this script (15 ABC-SMC sequences) with 10 cores
is approximately __*36 hours*__.
