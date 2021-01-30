% ECE302 -- Proj1
% Steven Lee & Jonathan Lam
% Prof. Keene
% 1/29/21

%%
clc; clear;
N = 100000;			% number of trials
D = 3;				% number of dice to take the sum of
T = 3;				% number of trials (in fun method)
S = 6;				% number of skills

%% 1a
q1a_exact = 1 / 6^D;

rolls = roll_dice([N D]);	% random simulation
q1a_est = sum(sum(rolls, 2) == 18) / N;

print_res('1a', q1a_exact, q1a_est);

%% 1b
q1b_exact = 1 - (1/216)^0*(215/216)^T;
rolls = roll_dice([N*T D]);
q1b_est = sum(sum(reshape(sum(rolls, 2) == 18, N, T), 2) >= 1) / N;

print_res('1b', q1b_exact, q1b_est);

%% 1c
q1c_exact = (q1b_exact)^S;
rolls = roll_dice([N*T*S D]);
q1c_est = sum(sum(reshape(sum(reshape(sum(rolls, 2) == 18, [], T), 2) >= 1, [], S), 2) == S) / N;

print_res('1c', q1c_exact, q1c_est);

%% 1d 
% P(all <= 9) - P(all <= 9 && none == 9)
q1d_exact = ((81/216)^T - ((81-25)/216)^T)^S;
 
rolls = roll_dice([N*T*S D]);
q1d_est = sum(sum(reshape(max(reshape(sum(rolls, 2), [], T), [], 2) == 9, ...
	[], S), 2) == S) / N;

print_res('1d', q1d_exact, q1d_est);

% helper function
function print_res(question, exact, est)
	pct_err = abs((est - exact) / exact * 100);
	fprintf('Question %s: Exact: %f; Estimate: %f; %% Err: %f%%\n', ...
		question, exact, est, pct_err);
end

% usage: roll_dice([x y])
function res = roll_dice(shape)
	res = randi(6, shape);
end