# GPLFM

A MATLAB implementation of the Gaussian process latent force model (GPLFM). 

This code is a supplement to [**J. Zou, E. Lourens, A. Cicirello (2022). "Virtual sensing of subsoil strain response in monopile-based offshore wind turbines via Gaussian process latent force models."**](https://arxiv.org/abs/2207.05901) Please cite this work if you use or refer to any parts of this repository. 


## Version history 

**Version 0.1.0**

Provides essential functions for implementing the GPLFM for joint input-state estimation using an MDOF cantilever system, cast into a modally reduced-order state-space model with acceleration-only measurements. 


## Documentation

The GPFLM is a Bayesian joint input-state estimator which infers the dynamic response of a structural system subject to unknown excitation. An augmented state-space model, which combines a physics-driven model of the system with a data-driven model of latent variables, characterizes the joint relationship between unknown inputs and states with a GP prior. Posterior inference of latent states is obtained via Gaussian process regression of measured states, performed sequentially using Kalman filtering and RTS smoothing.

**Getting started**

Run the file `run_GPLFM_Modal.m` in MATLAB (R2021a).


**Modules**

* `model` : functions which specify geometric and material properties (mass, stiffness, damping) of the structural system  

* `load` : functions which define the form, magnitude, and location of externally applied forces/inputs  

* `ss` : functions which assemble the state-space matrices of the structural model, including the modally reduced-order formulation  

* `rs` : functions which simulate dynamic response by numerical integration of EOM

* `kernel` : functions which assemble the state-space matrices of a GP expressed as a linear time-invariant stochastic differential equation based on a chosen GP ovariance kernel

* `hp` : functions which solve for hyperparameters of the GP covariance kernel  

* `id` : functions which perform joint input-state estimation by Kalman filtering/smoothing

* `plot` : functions for producing figures



**Example**

The provided numerical example implements the GPLFM to predict modal forces and response states (displacements, velocities, accelerations) of a 10-dof uniform cantilever structure with point masses $m$ and element stiffnesses $k$. The structure is expressed with a reduced-order model with $n_{mode} = 3$. Mass proportional damping with a modal damping ratio of $\xi$ is assumed at all degrees of freedom. 

An artificial input is defined as the combination of a sinusoidal force applied at the top (10th) level of the model and a Gaussian process with a Matern $\nu = 5/2$ covariance kernel with parameters $ \{ \alpha, l_s \}$, applied to all levels of the model. The synthetic system is constructed to loosely represent environmental and operational load conditions on a wind turbine support structure, where the sinusoidal force represents the harmonic frequency of vibration induced by the rotor motion of the turbine and the GP force represents ambient wind excitation distributed along the height of the tower. The input is non-Gaussian due to the contribution of the sinusoidal force.


![input](/figures/input.png)

<img src="https://github.com/joannajzou/GPLFM/tree/main/figures/input.png" width="200" />


The ground truth response of the cantilever structure to the artificial load is simulated using the Newmark average acceleration method. It is assumed that acceleration states at levels 5 and 10 are available as measurement data. Hyperparameters $ \{ \alpha, l_s \}$ of the GP covariance kernel are fit using the method in Section 2.3.1 of J. Zou et al. (2022), which minimizes the Hellinger distance between the empirical normal distribution fit to the measurement data and modeled GP prior corresponding to the hyperparameters. It is observed that the optimal hyperparameters from the tuning process lead to a strong match between the empirical and modeled distributions over the measurement data.

![hptuning](/figures/hptuning.png)
![hpfit](/figures/hpfit.png)


The GPLFM produces posterior estimates of response states, shown for level 5 of the structure, as well as of the first three modal forces. Since only $n_{mode} = 3$ modes are retained in the reduced-order model, accuracy in the frequency spectrum declines past the third natural frequency of the system. 

![responseestimation](/figures/responseestimation.png)
![modalforceestimation](/figures/modalforceestimation.png)




## Primary references

[1] [R. Nayek, S. Chakraborty, S. Narasimhan (2019). "A Gaussian process latent force model for joint input-state estimation in linear structural systems."](https://www.sciencedirect.com/science/article/abs/pii/S0888327019302286)

[2] [J. Hartikainen, S. Sarkka (2011). "Sequential inference for latent force models."](https://arxiv.org/abs/1202.3730)

[3] [J. Hartikainen, S. Sarkka (2010). "Kalman filtering and smoothing solutions to temporal Gaussian process regression models."](https://ieeexplore.ieee.org/document/5589113)





