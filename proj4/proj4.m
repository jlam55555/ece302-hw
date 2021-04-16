%%% ECE302 Project 4
% Steven Lee & Jonathan Lam
set(0,'defaultTextInterpreter','latex');
clc; clear; close all;

%% q1a
N = 1e4;

% prior probabilities
p0 = 0.8;

% choose variance
sigma = 0.5;

% choose constant value for A is a
a = 1;

% A is the signal, X is the gaussian additive noise
A = a*(rand(N,1) > 0.8);
X = sigma*randn(N, 1);
Y = X + A;

% MAP decision boundary
dec_boundary = (2*sigma^2*log(p0/(1-p0)) + a^2)/(2*a);

accuracy = mean((Y > dec_boundary) == A);

% TODO: write out formula for theoretical and compare results

%% q1b
etas = linspace(-10, 10, 1e3);
sigmas = logspace(-1, 1, 5);
P_F = zeros(length(etas), 1);
P_D = zeros(length(etas), 1);

% iterate over decision boundary
% and iterate over several SNRs
for j=1:length(sigmas)
    X = sigmas(j)*randn(N, 1);
    Y = X + A;

    for i=1:length(etas)
        % accuracy = mean((Y > etas(i)) == A);
        A_hat = Y > etas(i);

        P_D(i) = sum((A_hat == 1) & (A == 1)) / sum(A == 1);
        P_F(i) = sum((A_hat == 1) & (A == 0)) / sum(A == 0);
    end

    figure();
    plot(P_F, P_D);
    ylabel('$$P_D$$');
    xlabel('$$P_F$$');
    title(sprintf('Receiver Operating Characteristic $$\\sigma=%f$$', ...
        sigmas(j)));
end