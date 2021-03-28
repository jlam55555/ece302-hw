% ECE302 -- Project 4
% Steven Lee & Jonathan Lam
clc; clear; close all;

%% 1

% ML estimate function for alpha in rayleigh RV
est_fn_ray = @(samples) sqrt(mean(samples.^2)/2);

% values to loop over
Ns = logspace(1, 5, 5);
params = 1:3;

% dimensions: (value of N) x (value of param) x (exp/ray) x (bias/var/mse)
% results = zeros(length(Ns), length(params), 2, 3);

figure();
t1 = tiledlayout(3, length(params));

figure();
t2 = tiledlayout(3, length(params));

for param=params
    exp_res = zeros(length(Ns), 3);
    ray_res = zeros(length(Ns), 3);
    
    for Ni=1:length(Ns)
        for M=100
            N = Ns(Ni);
            
            [exp_bias, exp_var, exp_mse] = ...
                run_experiment(N, M, param, @exprnd, @mean);
            [ray_bias, ray_var, ray_mse] = ...
                run_experiment(N, M, param, @raylrnd, est_fn_ray);
            
            exp_res(Ni, :) = [exp_bias exp_var exp_mse];
            ray_res(Ni, :) = [ray_bias ray_var ray_mse];

            fprintf("N=%d param=%d: %f %f %f; %f %f %f\n", ...
                N, param, ...
                exp_bias, exp_var, exp_mse, ...
                ray_bias, ray_var, ray_mse);
        end
    end

    % for exponential
    nexttile(t1);
    loglog(Ns, abs(exp_res(:,1)));
    title('(Magnitude of) Bias');
    nexttile(t1);
    loglog(Ns, exp_res(:,2));
    title('Variance');
    nexttile(t1);
    loglog(Ns, exp_res(:,3));
    title('MSE');
    
    % for rayleigh
    nexttile(t2);
    loglog(Ns, abs(ray_res(:,1)));
    title('(Magnitude of) Bias');
    nexttile(t2);
    loglog(Ns, ray_res(:,2));
    title('Variance');
    nexttile(t2);
    loglog(Ns, ray_res(:,3));
    title('MSE');
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