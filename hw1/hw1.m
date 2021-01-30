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
F = 6;				% number of faces on the die

%% 1a
q1a_exact = 1 / 6^D;

rolls = roll_dice([N D], F);	% random simulation
q1a_est = sum(sum(rolls, 2) == 18) / N;

print_res('1a', q1a_exact, q1a_est);

%% 1b
q1b_exact = 1 - (1/216)^0*(215/216)^T;
rolls = roll_dice([N*T D], F);
q1b_est = sum(sum(reshape(sum(rolls, 2) == 18, N, T), 2) >= 1) / N;

print_res('1b', q1b_exact, q1b_est);

%% 1c
q1c_exact = (q1b_exact)^S;
rolls = roll_dice([N*T*S D], F);
q1c_est = sum(sum(reshape(sum(reshape(sum(rolls, 2) == 18, [], T), 2) >= 1, ...
	[], S), 2) == S) / N;

print_res('1c', q1c_exact, q1c_est);

%% 1d 
% P(all <= 9) - P(all <= 9 && none == 9)
q1d_exact = ((81/216)^T - ((81-25)/216)^T)^S;
 
rolls = roll_dice([N*T*S D], F);
q1d_est = sum(sum(reshape(max(reshape(sum(rolls, 2), [], T), [], 2) == 9, ...
	[], S), 2) == S) / N;

print_res('1d', q1d_exact, q1d_est);

%% 2a
D_T = 1;			% number of dice trolls use
F_T = 4;			% number of faces on troll die
D_K = 2;			% number of dice Keene uses
F_K = 2;			% number of faces on Keene die

% average hitpoints for trolls
q2a_exact = (1*1 + 2*1 + 3*1 + 4*1) / 4;
rolls = roll_dice([N, D_T], F_T);
q2a_est = mean(sum(rolls, 2));
print_res('2a', q2a_exact, q2a_est);

% average damage from fireball
q2a_exact = (2*1 + 3*2 + 4*1) / 4;
rolls = roll_dice([N, D_K], F_K);
q2a_est = mean(sum(rolls, 2));
print_res('2a', q2a_exact, q2a_est);

% TODO: bound on the probability?
% probability of dealing more than 3 damage with a fireball
q2a_exact = 1/4;
rolls = roll_dice([N, D_K], F_K);
q2a_est = sum(sum(rolls, 2) > 3) / N;
print_res('2a', q2a_exact, q2a_est);

%% 2b
q2b_exact = ones(1, 4) / 4;
rolls = roll_dice([N D_T], F_T);
q2b_est = histcounts(sum(rolls, 2), 'BinMethod', 'integers') / N;
fprintf('2b\n%d\n');
q2b_exact
q2b_est

q2b_exact = [0.25 0.5 0.25];
rolls = roll_dice([N D_K], F_K);
q2b_est = histcounts(sum(rolls, 2), 'BinMethod', 'integers') / N;
fprintf('2b\n');
q2b_exact
q2b_est

%% 2c 
% given the trolls, life, we willl evaluate whether or not the
% fireball can kill the trolls then we do that trial 6 times
%given troll health is 1 or 2, trolls always die, if health is
% 3 then 75% and if health is 4, then 25%
q2c_exact = 0.25*(1/2)^6+0.5*(3/4)^6+0.25;
rolls_troll = reshape(roll_dice([N*6 D_T], F_T), [] ,6);
rolls_fire = sum(roll_dice([N D_K], F_K),2);

diff = rolls_troll - rolls_fire;
maximum = max(diff, [], 2);
q2c_est = 1- sum(maximum > 0)/N;

print_res('2c', q2c_exact, q2c_est);

%% 2d
% scenarios for 5 dying, one surviving:
% note that it is only possible that FB=2 or 3 and HP=3 or 4
% define the following events:
% a := HP=4, FB=3, all other HP<=3; P(a) = (1/4)(1/2)(3/4)^5
% b := HP=4, FB=2, all other HP<=2: P(b) = (1/4)(1/4)(1/2)^5
% c := HP=3, FB=2, all other HP<=2: P(c) = (1/4)(1/4)(1/2)^5 = P(b)
% P(d) := P(HP=4 ^ 5 dying, one surviving) = P(a) + P(b)
% P(e) := P(HP=3 ^ 5 dying, one surviving) = P(c)
% P(f) := P(5 dying, one surviving) = P(a ^ b ^ c) = P(a) + P(b) + P(c)
% E[HP | 5 dying, 1 surviving] = 4*P(d)/P(f) + 3*P(e)/P(f) (+ 2*0 + 1*0)
Pa = 1/4 * 1/2 * (3/4)^5;
Pb = 1/4 * 1/4 * (1/2)^5;
Pc = Pb;
Pd = Pa + Pb;
Pe = Pc;
Pf = Pa + Pb + Pc;

q2d_exact = 4*Pd/Pf + 3*Pe/Pf;

rolls_T = roll_dice([N*6 D_T], F_T);
rolls_K = roll_dice([N D_K], F_K);

% in simulation, need to select the cases where exactly 5 die
troll_hp = reshape(sum(rolls_T, 2), [], 6);
fire_dmg = sum(rolls_K, 2);
matches = troll_hp - fire_dmg;
matches_where_1_troll_wins = sum(matches > 0, 2) == 1;
matches = matches(matches_where_1_troll_wins, :);
troll_hp = troll_hp(matches_where_1_troll_wins, :);
q2d_est = mean(troll_hp(matches > 0));

print_res('2d', q2d_exact, q2d_est);

% helper function
function print_res(question, exact, est)
	pct_err = abs((est - exact) / exact * 100);
	fprintf('Question %s: Exact: %f; Estimate: %f; %% Err: %f%%\n', ...
		question, exact, est, pct_err);
end

% usage: roll_dice([x y])
function res = roll_dice(shape, faces)
	res = randi(faces, shape);
end