# Black-Sea_Bayes
This respiratory contains the code to back-track plastic particles from their sampling location to their more probable sources by using Bayesian inference, coupled to Lagrangian simulations of virtual particles.

##Simulations

The simulations are carried out by using the Ocean Parcels framework (version 2.0) which is fully describedin Lange and van Sebille (2017) and Delandmeter and van Sebille (2019).  The Ocean Parcels framework is freely available at <http://oceanparcels.org>.

To run the simulations use Simulation/Advection.py. However, the grid file needs to be created beforehand, by using Simulation/ParticleGrid.py.

##Entropy
The entropy is used in order to understand over which timescales the back-tracking of plastic particles can be carried over.

##Bayesian inference

Use Bayes_BlackSea.py
