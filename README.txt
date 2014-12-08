This project is intended to solve a set of reactions that controls daytime variations of ozone in the troposphere using forward Euler, backward Euler and Multistep Implicit-Explicit numerical solutions. Comparison of model output with observational data from TORCH 2 in WAO (UEA) is included.
The set of reactions is the following:
O -> O3 (k1)
NO2 + hv -> NO + O (J)
NO + O3 -> NO2 + O2 (k3)
This is so called null cycle, with no net formation of O3.