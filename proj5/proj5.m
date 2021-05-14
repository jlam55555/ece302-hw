%% ECE302 Project 5
%% Stevel Lee & Jonathan Lam

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
% 

% We have tried to replicate the system using the autocorrelation
% (Rrr) as well as cross-correlation (Rsr)

% So we used the xcorr function is matlab 
% (xcorr(conv(s, c, 'same') + d, s, 4, 'biased')) 
% and the results seems right, however, the mse is off the mark. 
% The accuracy look to be on mark, which even though the MSE does not





N = 10;
sig2 = 1;

% R_rr = [1.2+sig2 .28 .4 0 0 0 0 0 0 0];
R_rr = [1.2 .28 .4 0 0 0 0 0 0 0];
R = R_rr(abs((1:N) - (1:N).') + 1)
R_sr = [1 .2 .4 0 0 0 0 0 0 0].'

h = [0; 0; 0; 0; 0; 0; 0; 0; 0; (R \ R_sr(1:N))]

% h = [.7754; -.1410; -.305; .1175]

%% Part 2
M = 1e6;
s = 2*randi(2, [M 1]) - 3;
d = sqrt(sig2) * randn([M 1]);

c = [0 0 1 .2 .4];

r = conv(s, c, 'same'); % + d;
s_hat = conv(r, h, 'same');

mse = mean((s-s_hat).^2)


s(1:10)
d(1:10)

rnge = 100:200;
[s(rnge) sign(s_hat(rnge)) (s(rnge) == sign(s_hat(rnge))) s_hat(rnge)]
[mean((s(rnge) == sign(s_hat(rnge)))) mse]