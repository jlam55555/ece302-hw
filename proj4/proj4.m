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

%% q1b&c
etas = linspace(-10, 10, 1e3);
sigmas = logspace(-1, 1, 5);
P_F = zeros(length(etas), 1);
P_D = zeros(length(etas), 1);

% iterate over decision boundary
% and iterate over several SNRs
for j=1:length(sigmas)
    X = sigmas(j)*randn(N, 1);
    Y = X + A;

    % MAP rule: minimizing probability of error
    map_boundary = (2*sigma^2*log(p0/(1-p0)) + a^2)/(2*a);

    % q1c
    % 8.8 in detection theory pdf
    % (C01 - C11)*P1*f(y|H1) = (C10 - C00)*P0*f(y|H0)
    % Now set C11=C00=0 (as before), but set C01=10*C10 (rather than C01=C10)
    %
    % basically changes the coefficients from (P1,P0) to (10*P1, P0)
    % see new factor of 10
    q1c_boundary = (2*sigma^2*log(p0/(10*(1-p0))) + a^2)/(2*a);

    for i=1:length(etas)
        A_hat = Y > etas(i);

        P_D(i) = sum((A_hat == 1) & (A == 1)) / sum(A == 1);
        P_F(i) = sum((A_hat == 1) & (A == 0)) / sum(A == 0);

        % if this is the closest point to the decision boundary for part c
        if abs(etas(i) - q1c_boundary) < (etas(2) - etas(1)) / 2
            q1c_i = i;
        end

        if abs(etas(i) - map_boundary) < (etas(2) - etas(1)) / 2
            map_i = i;
        end
    end

    figure();
    hold('on');
    plot(P_F, P_D);
    plot(P_F(map_i), P_D(map_i), 'g*');
    plot(P_F(q1c_i), P_D(q1c_i), 'r*');
    hold('off');
    ylabel('$$P_D$$');
    xlabel('$$P_F$$');
    title(sprintf('Receiver Operating Characteristic $$\\sigma=%f$$', ...
        sigmas(j)));
    legend(["ROC", sprintf('MAP boundary (\\eta=%f)', map_boundary), ...
        sprintf('Q1c modified cost boundary (\\eta=%f)', q1c_boundary)]);
end