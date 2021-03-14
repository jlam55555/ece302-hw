% ECE302 Proj 2
% (L)MMSE estimators
% Jonathan Lam & Steven Lee
clc; clear; close all;
set(0, 'defaultTextInterpreter', 'latex');

%% 1a
N = 1e6;

Y = 2*rand(N, 1) - 1;
W = 4*rand(N, 1) - 2;

X = Y + W;

% piecewise estimators
X1 = X < -1;
X2 = (X >= -1) & (X <= 1);
X3 = X > 1;
Y1 = (1 + X(X1)) / 2;
Y2 = zeros(size(X(X2)));
Y3 = (-1 + X(X3)) / 2;

X = [X(X1); X(X2); X(X3)];
Y = [Y(X1); Y(X2); Y(X3)];
Y_est = [Y1; Y2; Y3];
Y_err = (Y - Y_est).^2;

% for plotting empirical piecewise MSE as a function of X
Y_X_err = sortrows([X Y_err]);
Y_X_err(:,2) = movmean(Y_X_err(:,2), N/1e2);

figure('Position', [0 0 1000 1000]);
tiledlayout(1, 1, 'TileSpacing', 'tight');
hold on;

% plot all X vs error data
scatter(X, Y_err);

% plot piecewise MSE
x = -3:1e-2:-1;
plot(x, (3+x).^2/12, 'k');
x = -1:1e-2:1;
plot(x, 1/3*ones(size(x)), 'k', 'HandleVisibility', 'off');
x = 1:1e-2:3;
plot(x, (3-x).^2/12, 'k', 'HandleVisibility', 'off');

% plot empirical piecewise MSE
plot(Y_X_err(:,1), Y_X_err(:,2), 'r');

legend({'Squared Error', 'Theoretical MSE', 'Empirical MSE (moving average)'});
xlabel('Observed $$X$$');
ylabel('$$(\hat{Y}(X)-Y)^2$$');
title('Error in estimate vs. observed $$X$$');

% calculate empirical MSE
MSE_emp = mean(Y_err);

%% 1b
Y_lin = X/5;
LMSE_emp = mean((Y_lin - Y).^2);

%table making
display(table([1/4; 4/15], [MSE_emp; LMSE_emp], ...
	'VariableNames', {'Theoretical', 'Empirical'}, ...
	'RowNames', {'MMSE', 'LMMSE'}))

%% 2
clear;
N = 1e4;
mu_y = 1;
sig_ys = 1:10;
sig_rs = 1:10;
Rs = 1:10;		% R is the number of observations

% calculate all da results
res = [];
row_labels = {};
i = 1;
for sig_y_ind = 1:length(sig_ys)
for sig_r_ind = 1:length(sig_rs)
for R_ind = 1:length(Rs)
	sig_y = sig_ys(sig_y_ind);
	sig_r = sig_rs(sig_r_ind);
	R = Rs(R_ind);
	[emp, theo] = calc_llmse(mu_y, sig_y, sig_r, N, R);
	res(R_ind, sig_r_ind, sig_y_ind, :) = [emp, theo];
	row_labels{i} = sprintf('sig_y=%.2f, sig_r=%.2f, R=%d', sig_y, sig_r, R);
	i = i + 1;
end
end
end

% plotting slices of the result vector
figure('Position', [0 0 1200 400]);
tiledlayout(1, 3, 'TileSpacing', 'tight');
% plot MSE vs R
nexttile();
mvr = reshape(res(:, 5, 5, :), [], 2);
plot(Rs, mvr(:,1), Rs, mvr(:,2));
xlabel('Number of observations');
ylabel('LMMSE');
title('LMMSE vs. Number of observations ($$\sigma_y=\sigma_r=5$$)');
legend({'Empirical LMMSE', 'Theoretical LMMSE'});
% plot MSE vs sig_y
nexttile();
mvsy = reshape(res(5, :, 5, :), [], 2);
plot(sig_ys, mvsy(:,1), sig_ys, mvsy(:,2));
xlabel('$$\sigma_y$$');
ylabel('LMMSE');
title('LMMSE vs. $$\sigma_y$$ ($$R=\sigma_r=5$$)');
legend({'Empirical LMMSE', 'Theoretical LMMSE'});
% plot MSE vs sig_r
nexttile();
mvsr = reshape(res(5, 5, :, :), [], 2);
plot(sig_rs, mvsr(:,1), sig_rs, mvsr(:,2));
xlabel('$$\sigma_r$$');
ylabel('LMMSE');
title('LMMSE vs. $$\sigma_r$$ ($$R=\sigma_y=5$$)');
legend({'Empirical LMMSE', 'Theoretical LMMSE'});

% displaying results in a table
res = reshape(res, [], 2);
display(table(res(:,1), res(:,2), ...
	'VariableNames', {'Empirical', 'Theoretical'}, 'RowNames', row_labels))

% function to help with calculations at different values of the parameters
function [emp, theo] = calc_llmse(mu_y, sig_y, sig_r, N, R)
	Y_gaus = mu_y + sig_y * randn(N, 1);
	R_gaus = sig_r * randn(N, R);
	X = Y_gaus + R_gaus;

	% solve normal equation C_XX*a=C_XY (8.70)
	C_XX = sig_y^2 * ones(R) + sig_r^2 * eye(R);
	C_XY = sig_y^2 * ones(R, 1);
	a = C_XX^-1 * C_XY;

	% 8.74: comes from requirement that Yhat_lin is unbiased
	% (same as 8.60 assuming that mu_x=mu_y since mu_r=0)
	a0 = mu_y * (1 - sum(a));

	Y_gaus_est = a0 + sum(a.' .* X, 2);

	% empirical MSE
	emp = mean((Y_gaus_est - Y_gaus).^2);

	% 8.73: theoretical associated MSE
	theo = sig_y^2 - C_XY.' * C_XX^-1 * C_XY;
end