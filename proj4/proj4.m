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

        P_D(i) = sum((A_hat == a) & (A == a)) / sum(A == a);
        P_F(i) = sum((A_hat == a) & (A == 0)) / sum(A == 0);

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

%% q1d
% use same sigma as in part a
X = sigma*randn(N, 1);
Y = X + A;

% iterate over values of the prior probability
p0s = linspace(0, 1, 1e2);
costs = zeros(length(p0s), 1);

for i=1:length(p0s)
    % find best decision boundary
    q1c_boundary = (2*sigma^2*log(p0s(i)/(10*(1-p0s(i)))) + a^2)/(2*a);
    A_hat = Y > q1c_boundary;

    % get cost
    costs(i) = 10 * mean((A_hat == 0) & (A == a)) ...   % false negative
        + mean((A_hat == a) & (A == 0));                % false positive
end

figure();
plot(p0s, costs);
ylabel('Mean cost');
xlabel('$$P_0$$ (prior probability of target not present)');
title('Cost vs. prior probability of target not present');

%% q1e
% TODO: later

%% q2
clc; clear; close all;
load('Iris');

N = size(features, 1);          % numbere of samples
C = length(unique(labels));     % number of classes
K = size(features, 2);          % number of features

% split into train/test sets
% train-test split 50/50
is_train        = rand(N, 1) < 0.5;
train_features  = features(is_train, :);
train_labels    = labels(is_train, :);
test_features   = features(is_train == 0, :);
test_labels     = labels(is_train == 0, :);

% store results of evaluating model on test dataset
results = zeros(length(test_labels), C);

% calculate class priors (on train dataset)
priors = histcounts(test_labels) / length(test_labels);

for i=1:C
    indices = train_labels == i;
    class_features = train_features(indices, :);
    
    % "train": find MAP parameters (class-conditional density)
    mus = mean(class_features);
    covs = cov(class_features);
    
    % evaluate on test set
    results(:,i) = mvnpdf(test_features, mus, covs) * priors(i);
end

[mx, est] = max(results, [], 2);
accuracy = mean(est == test_labels)