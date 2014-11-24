This program is intended to solve a set of reactions that affects ozone production and destruction in the tropo- or stratosphere at daytime using Forward Euler, Backward Euler and Multistep Implicit-Explicit numerical schemes.
The initial set of reactions is the following:
NO + O3 -> NO2 + O2 (k1)
NO2 + hv -> NO + O (J)
O -> O3 (k3)
This is a null cycle with no net formation of O3.