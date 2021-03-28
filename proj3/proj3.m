% ECE302 -- Project 4
% Steven Lee & Jonathan Lam
clc; clear; close all;

%% 1

% ML estimate function for alpha in rayleigh RV
est_fn_ray = @(samples) sqrt(mean(samples.^2)/2);

for N=logspace(1, 5, 5)
for M=100
for param=1
	[exp_bias, exp_var, exp_mse] = ...
        run_experiment(N, M, param, @exprnd, @mean);
	[ray_bias, ray_var, ray_mse] = ...
        run_experiment(N, M, param, @raylrnd, est_fn_ray);

	fprintf("N=%d param=%d: %f %f %f; %f %f %f\n", ...
		N, param, ...
		exp_bias, exp_var, exp_mse, ...
		ray_bias, ray_var, ray_mse);
end
end
end


% N = number of samples
% M = number of trials
% param = actual parameter
% randfn = function to generate random samples
% estfn = calculate ML estimate of parameter given samples
function [bias, variance, mse] = run_experiment(N, M, param, randfn, estfn)
    % generate M samples of N values sampled from the distribution
	samples = randfn(param, N, M);
    
    % generate ML estimate of variable
	est = estfn(samples);

	% calculate bias, variance, MSE (eq. 8.18) of estimator
	bias = mean(est) - param;
	variance = var(est);
	mse = variance + bias^2;
end