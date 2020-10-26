# Black-Sea_Bayes
This respiratory contains the code to back-track plastic particles from their sampling location to their more probable sources by using Bayesian inference, coupled to Lagrangian simulations of virtual particles.

##Simulations
The simulations are carried out by using the Ocean Parcels framework (version 2.0) which is fully describedin Lange and van Sebille (2017) (article: <https://arxiv.org/abs/1707.05163>) and Delandmeter and van Sebille (2019) (article: <https://dspace.library.uu.nl/handle/1874/384225>). The Ocean Parcels framework is freely available at <http://oceanparcels.org>.

To run the simulations use Simulation/Advection.py. However, the grid file needs to be created beforehand, by using Simulation/ParticleGrid.py.

##Entropy
The entropy is used in order to understand over which timescales the back-tracking of plastic particles can be carried over as suggested by Wichmann et al. (2019) (Link to the article : <https://iopscience.iop.org/article/10.1088/2515-7620/ab4e77/meta>).

##Bayesian inference
Use Bayes/Bayes_BlackSea.py. Beforehand, the density of particles after x-time of simulation must be calculated by using Bayes/Density.py.

