# Bayesian-Reconciliation-of-the-Return-Predictability.

This repository contains code for: *"Bayesian Reconciliation of the Return Predictability"* (joint work with [Sylvia Frühwirth-Schnatter](https://statmath.wu.ac.at/~fruehwirth/) and [Leopold Sögner](https://myservice.ihs.ac.at/soegner/)).

We share the code to obtain simulated data and financial data. To receive similar results please change the number of **MC_draws**, **burnin** and **thining** as it's specified in the paper. We further provide code to replicate Figures and Tables in the paper.

**Code**:
1) start_IN.R - the main file that loads all required packages and functions (from functions_IN.R). Please run first.
2) functions_IN.R - the file contains the necessary functions to replicate the results.
3) sim_data_IN.R - the file contains the code to replicate the simulated data and to create the raw results.
4) real_data_IN.R - the file contains the code to construct the empirical data set, and perform Bayesian parameter estimation and testing.
5) oos_analysis_IN.R - the file contains the code to perform the OOS analysis.
6) figures_IN.Rmd - the file contains the code to replicate figures and tables provided in the paper.
7) figures_IN.html - the output of the figures_IN.Rmd

**Data**:
It should be downloaded and added to the main folder
https://www.dropbox.com/scl/fo/uw91aktmh1p7ond33uaqk/h?rlkey=vly8j10pxmr0lzbrmeex4n37g&dl=0
1) real_data - the folder that contains the raw real data and the results obtained for the real-world data.
2) sim_data - the folder containing the simulated data results.

**Supplementary** - the folder that contains the results for Table 9 where we used the EViews package.

Please cite the corresponding paper:
```
@article{fsk_bayesian,
  title={Bayesian Reconciliation of Return Predictability},
  author={Koval, Borys and Fr{\"u}hwirth-Schnatter, Sylvia and S{\"o}gner, Leopold},
  journal={arXiv preprint arXiv:2212.02239},
  year={2022}
}
```
or visit [SSRN page](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4288973). In case of any questions please contact Borys [borys.koval@wu.ac.at](mailto:borys.koval@wu.ac.at)
