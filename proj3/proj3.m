% ECE302 -- Project 4
% Steven Lee & Jonathan Lam
clc; clear; close all;
set(0, 'defaultTextInterpreter', 'latex');

% ML estimate function for alpha in rayleigh RV
est_fn_ray = @(samples) sqrt(mean(samples.^2)/2);

%% Q1

% values to loop over
Ns = logspace(1, 5, 5);
params = 1:4;

% one plot for exponential, one for rayleigh
figure('Position', [0 0 750 1000]);
t1 = tiledlayout(length(params), 3);
sgtitle("Exponential R.V.");

figure('Position', [0 0 750 1000]);
t2 = tiledlayout(length(params), 3);
sgtitle("Rayleigh R.V.");

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
    semilogx(Ns, exp_res(:,1));
    title(sprintf('Bias ($$\\mu=%d$$)', param));
    xlabel("Samples");
    ylabel("$$E[\hat\mu-\mu]$$");
    grid on;
    
    nexttile(t1);
    loglog(Ns, exp_res(:,2));
    title('Variance');
    xlabel("Samples");
    ylabel("$$E[(\hat\mu-\bar{\hat\mu})^2]$$");
    grid on;
    
    nexttile(t1);
    loglog(Ns, exp_res(:,3));
    title('MSE');
    xlabel("Samples");
    ylabel("$$E[(\hat\mu-\mu)^2]$$");
    grid on;
    
    % for rayleigh
    nexttile(t2);
    semilogx(Ns, ray_res(:,1));
    title(sprintf('Bias ($$\\alpha=%d$$)', param));
    xlabel("Samples");
    ylabel("$$E[\hat\alpha-\alpha]$$");
    grid on;
    
    nexttile(t2);
    loglog(Ns, ray_res(:,2));
    title('Variance');
    xlabel("Samples");
    ylabel("$$E[(\hat\alpha-\bar{\hat\alpha})^2]$$");
    grid on;
    
    nexttile(t2);
    loglog(Ns, ray_res(:,3));
    title('MSE');
    xlabel("Samples");
    ylabel("$$E[(\hat\alpha-\alpha)^2]$$");
    grid on;
end

% Explanation of figures:
% - The first figure is for the exponential R.V., the other for Rayleigh.
% - Each column represents one value of the parameter (mu or lambda).
% - For each value of the parameter, multiple (10) trials of different
%   sample sizes (N={10,100,1000,10000,100000}) were used. The bias,
%   variance, and MSE of the estimator was taken from these.
% - Each plot plots N (number of samples) on the x-axis vs. the specified
%   metric.


%% Q2
load('data');
data = data.';      % need a column vector

% range of values of data
x = linspace(min(data), max(data), 1000);

mu_est = mean(data);
alpha_est = est_fn_ray(data);

lambda_est = 1/mu_est;
exp_pdf = @(x) lambda_est * exp(-lambda_est * x);

alpha2_est = alpha_est^2;
ray_pdf = @(x) x/alpha2_est .* exp(-x.^2/(2*alpha2_est));

% see which distribution generates the higher likelihood
% (use sum of log-likelihoods, also can use product of likelihoods)
exp_likelihood = sum(log(exp_pdf(data)));
ray_likelihood = sum(log(ray_pdf(data)));

% rayleigh has higher likelihood, so it is the more likely distribution
fprintf("Log-likelihood of exponential distribution: %f\n" ...
    + "Log-likelihood of rayleigh distribution: %f\n", ...
    exp_likelihood, ray_likelihood);

%% show histogram to visually see which distribution fits better
figure();
hold on;
histogram(data, 'Normalization', 'pdf');
plot(x, exp_pdf(x));
plot(x, ray_pdf(x));
legend(["Sample data (normalized to PDF)", ...
    sprintf("Exponential PDF (mu=%f)", mu_est), ...
    sprintf("Rayleigh PDF (alpha=%f)", alpha_est)]);
title("Data vs. Exponential and Rayleigh ML-Estimated Distributions");
ylabel("PDF");
xlabel("Values");

% we see that the histogram closely matches the Rayleigh distribution,
% so it most likely is drawn from this distribution


%% helper function to run experiment

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
    
    bias = mean(est) - param;
    variance = var(est);
    mse = variance + bias^2;
end