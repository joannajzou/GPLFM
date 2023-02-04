# GPLFM

A MATLAB implementation of the Gaussian process latent force model (GPLFM). 

The GPFLM a Bayesian joint input-state estimator which infers the dynamic response of a structural system subject to unknown excitation by endowing latent variables with a GP prior. Posterior inference of latent states is obtained via Gaussian process regression of observed states, performed sequentially using Kalman filtering and RTS smoothing with an augmented state-space model.

This code is a supplement to [J. Zou, E. Lourens, A. Cicirello (2022). "Virtual sensing of subsoil strain response in monopile-based offshore wind turbines via Gaussian process latent force models."](https://arxiv.org/abs/2207.05901) Please cite this work if you use or refer to any parts of this repository. 


## Version history 
**Version 0.1.0:** Provides essential functions for implementing the GPLFM for joint input-state estimation using an MDOF cantilever system, cast into a modally reduced-order state-space model with acceleration-only measurements. 


## Documentation

run_GPLFM_Modal.m : main MATLAB script

**Modules**
    model : functions which specify geometric and material properties (mass, stiffness, damping) of the structural system
    load : functions which define the form, magnitude, and location of externally applied forces/inputs
    ss : functions which assemble the state-space matrices of the structural model, including the modally reduced-order formulation
    rs : functions which simulate dynamic response (e. g. Newmark average acceleration method)
    kernel: functions which assemble the state-space matrices of a GP expressed as a linear time-invariant stochastic differential equation (SDE), based on a chosen covariance kernel function
    hp : functions which solve for hyperparameters of the GP covariance kernel
    id : functions which perform joint input-state estimation
    plot : functions for plotting results



## Primary references

[1] [R. Nayek, S. Chakraborty, S. Narasimhan (2019). "A Gaussian process latent force model for joint input-state estimation in linear structural systems."](https://www.sciencedirect.com/science/article/abs/pii/S0888327019302286)

[2] [J. Hartikainen, S. Sarkka (2011). "Sequential inference for latent force models."](https://arxiv.org/abs/1202.3730)

[3] [J. Hartikainen, S. Sarkka (2010). "Kalman filtering and smoothing solutions to temporal Gaussian process regression models."](https://ieeexplore.ieee.org/document/5589113)





