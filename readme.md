Examples for "How Kalman Filters Work"
======================================

This repository contains the MATLAB files used for the "How Kalman Filters
Work" article available at:

http://www.anuncommonlab.com/articles/how-kalman-filters-work/

including four filters, demonstration files, and several utilities.

Filters:

bf   bootstrap filter
ekf  extended Kalman filter
lkf  linear Kalman filter
ukf  sigma-point filter

Demos:

ekf_demo          demonstration of an extended Kalman filter
kalman_gain_demo  demonstration of the correction of the covariance via the
                  Kalman gain
lkf_demo          demonstration of a linear Kalman filter with comparisons
                  to the extended Kalman filter
sigma_point_dmeo  demonstration of a sigma-point (unscented) Kalman filter
particle_demo     demonstration of a bootstrap filter

Utilities:

covdraw  utility to create random draws from a given covariance matrix
ellipse  creates points of an ellipse corresponding to a covariance matrix
mndpdf   utility to calculate probability density for a multivariate normal
         distribution
rk4step  numerically integrates a function over one time step using the
         common Runge-Kutta fourth order method

Copyright 2016 An Uncommon Lab
