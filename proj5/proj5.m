%% ECE302 Project 5
%% Steven Lee & Jonathan Lam

clc;
clear;
close all;

%% Part 1

%                        d[n]
%                         |
%         +------+        v        +------+
% s[n] -> | c[n] | ----> (+) ----> | h[n] | -> shat[n]
%         +------+  x[n]      r[n] +------+
%
% Calculations: (note: * means convolution)
%
% R_ss[n] = E[s[m]s[m+n]] = { 1   if n=0 } = delta[n]
%                           { 0   else   }
%
% Use eq. 11.8 and fact that s, d are uncorrelated:
% R_rs[n] = R_xs[n] = R_ss[n] * c[n] = c[n]
%
% Use fact that s, d are uncorrelated:
% R_rr[n] = R_xx[n] + R_dd[n]
% 
% Autocorrelation of gaussian white noise:
% R_dd[n] = sigma^2 * delta[n]
% (we set sigma^2=1)
%
% based on Wiener-Khinchin theorem:
% R_xx[n] = R_ss[n] * R_cc[n]
% where R_cc[n] = autocorrelation of c[n] = c*c

% We have tried to replicate the autocorrelation calculation
% R_xy[n] using the xcorr function xcorr(x, y, N, 'biased'))
% and the autocorrelations seems correct.

%% Part 2

% arbitrary choose a value for sigma^2
sig2 = 0.5;

% M is number of samples
M = 1e6;

% create random vectors
s = 2*randi(2, [M 1]) - 3;
d = sqrt(sig2) * randn([M 1]);

% need to pad c so that convolution with
% type "same" works correctly
c = [0 0 1 .2 .4];

% precompute r -- doesn't depend on h
r = conv(s, c, 'same') + d;

% N is the length of the filter h
for N = [4 6 10]
	% set up normal equations
	R_rr = zeros([N, 1]);
	R_rr(1:3) = [1.2+sig2 .28 .4];
	
	R_sr = zeros([N, 1]);
	R_sr(1:3) = [1 .2 .4].';

	R = R_rr(abs((1:N) - (1:N).') + 1);

	% solve normal equation
	% Rh = R_sr => h = inv(R)*R_sr
	h = R \ R_sr(1:N);

	% need to pad h so that it's correctly centered
	% so that conv w/ "same" padding works correctly
	h = [zeros([N-1 1]); h];

	% calculate estimate with our filter
	s_hat = conv(r, h, 'same');

	% calculate and print mse
	mse = mean((s-s_hat).^2);
	acc = mean(sign(s) == sign(s_hat));
	fprintf('N=%d: MSE=%f; accuracy=%f\n', N, mse, acc);
end

% For some reason, increasing N does not improve MSE. However, the
% MSE is fairly low, and the accuracy (determined by the relative
% frequency that s_hat predicts the correct sign of s) is definitely
% above random guessing (0.5). This may be because the coefficients
% of h[n] for n>4 are relatively small, so increasing N does not
% improve performance. However, changing sigma^2 affects performance
% as expected.