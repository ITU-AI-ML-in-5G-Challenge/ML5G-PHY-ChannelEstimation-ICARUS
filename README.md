# ML5G-PHY-ChannelEstimation-ICARUS
This repository contains the source codes and the description of the proposed solution to the challenge [ML5H-PHY [channel estimation]](https://research.ece.ncsu.edu/ai5gchallenge/) by the team ICARUS. The members of ICARUS are Özlem Tuğfe Demir, Cenk M. Yetis, Emil Björnson, and Pontus Giselsson.

Sparse Bayesian Learning for Site-Specific Hybrid MIMO Channel Estimation
==================


The package contains the description file for the proposed sparse Bayesian learning method and source codes based on Matlab.


## Abstract 

To tackle this challenging problem, we propose
mainly a sparse Bayesian learning algorithm to exploit the
sparsity of the channel. We utilize the pattern-coupling concept
to model possible block sparsity patterns among the consecutive
angle-of-arrivals (AoAs) and angle-of-departures (AoDs). As a
first step, we obtain the time-domain channels from the provided
training dataset by inverse FFT and remove the negligibly small
taps. Then, we apply the algorithm to the true time-domain
channels to obtain the sparse representations. Using joint angular
distribution that is learned from training data, we refine the
grids and pattern-coupling relations in testing stage in an aim to
improve the channel estimation quality.


## Content 

The description of the system model and the details of the proposed solution are provided in the file "description_ICARUS.pdf". The code "learning_angular_distribution.m" runs the algorithm to estimate the sparse vectors from the training data by loading the training data that is provided in https://research.ece.ncsu.edu/ai5gchallenge/. The same code saves the estimated sparse vectors in ".mat" files, which can be found under the folder "estimated_sparse_vectors_from_training_data" in https://1drv.ms/u/s!AsV9iFv-q8_lv0Oasnwa-a-OOcnJ?e=Ev8Zs8. Next, the code in "analyze_training_data_final.m" uses those estimates to apply the grid construction and saves the learned parameters in the file "learned_parameters.mat". Finally, the code in "algorithm_on_test_data_final.m" loads those parameters and run the algorithm on the testing data that is provided in https://research.ece.ncsu.edu/ai5gchallenge/. The estimated channels can be accessed from the ".mat" files under the folder "estimated_channels_nine_mat_files" in https://1drv.ms/u/s!AsV9iFv-q8_lv0Oasnwa-a-OOcnJ?e=Ev8Zs8.



## Acknowledgements

The work of all team members was supported by the Wallenberg AI, Autonomous Systems and Software Program (WASP) funded by the Knut and Alice Wallenberg Foundation. 
