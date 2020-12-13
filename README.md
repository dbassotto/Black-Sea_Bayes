# Black-Sea_Bayes
This respiratory contains the code to back-track plastic particles from their sampling location to their more probable sources by using Bayesian inference, coupled to Lagrangian simulations of virtual particles.

## Simulations
The simulations are carried out by using the Ocean Parcels framework (version 2.0) which is fully describedin Lange and van Sebille (2017) (article: <https://arxiv.org/abs/1707.05163>) and Delandmeter and van Sebille (2019) (article: <https://dspace.library.uu.nl/handle/1874/384225>). The Ocean Parcels framework is freely available at <http://oceanparcels.org>.

To run the simulations use Advection.py. However, the grid file needs to be created beforehand, by using ParticleGrid.py.

## Entropy
The Shannon entropy is useful to understand over which timescales the back-tracking of plastic particles can be done, as suggested by Wichmann et al. (2019) (Link to the article : <https://iopscience.iop.org/article/10.1088/2515-7620/ab4e77/meta>).

Use BS_Entropy_Shannon.py.

## Bayesian inference
Use BS_Bayes.py. Beforehand, the density of particles after x-time of simulation must be calculated by using BS_Density.py.

