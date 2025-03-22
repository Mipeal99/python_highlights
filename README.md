TO DO, MAYBE: translate most of the comments in the programs to English
Here is the gist of what these programs are:

-CAHA_reduction_script.ipynb: this is an astronomical data reduction project. It takes calibration data and makes bias, dark and flats and applies it to science data of NGC2392, taken from the CAFOS telescope in Calar Alto.

-prestacking.py: second part of the data reduction project. Does some fixing on the reduced data to prepare them for exposure stacking.

-adaptive_integration.py: implementation of Simpson's method for numerical integration with an adaptive step in the integration interval based on the behaviour of the integrand.

-hydrogenatom_numerov.py: solves the Schrödinger equation for the hydrogen atom using Numerov's method and calculates some of the wave functions.

-lithiumatom_numerov.py: the same as hydrogenatom_numerov.py, but adapted for a lithium atom with a hydrogen-like approximation.

-idealgas_simulation.py: graphical simulation of a 2D ideal gas with molecules moving inside a box. Has temperature and pressure calculations at the end, as well as adjustable initial conditions.

-masterthesis_numericalsolutions.py: solves the equations of motion of some models for an axion-like-particle interacting with a photon, calculates the photon survival probability for some initial conditions inspired on the Perseus galaxy cluster, and compares the numerical solutions with the analytical solution of the model that is most commonly used in the literature (what I call "the G_0 model" in the program). There are also many commented sections that are mostly additional plots for the thesis and the presentation.

-masterthesis_constraints.py: contains most of the code of the previous program but it also generates a simulation of experimental data from the Perseus galaxy cluster (based on MAGIC results), modifies the simulation to account for the hypothetical effects of axion-like particles for a givien parameter set and performs a chi squared test to see if the chosen parameters are compatible with the error bars. The rejected points of the parameter space can be later summarized on an exclusion plot.

-maxwellsequations_1.py: solves Maxwell's equations in two dimensions and graphs a wave propagating in a box.

-maxwellsequations_2.py: solves Maxwell's equations in two dimensions and simulates two different mediums (allowing for different conductivities) separated by a double slit. The interference pattern is graphed in real time as the electromagnetic wave propagates.
