% clc; clear; close all;

%% 1a
N = 1e6;

Y = 2*rand(N, 1) - 1;
W = 4*rand(N, 1) - 2;

X = Y + W;

X1 = X < -1;
X2 = (X >= -1) & (X <= 1);
X3 = X > 1;

Y1 = (1 + X(X1)) / 2;
Y2 = zeros(size(X(X2)));
Y3 = (-1 + X(X3)) / 2;

YYhat = [Y(X1) Y1; Y(X2) Y2; Y(X3) Y3];

Xs = [X(X1); X(X2); X(X3)];
Ys = (YYhat(:,1) - YYhat(:,2)).^2;

alldata = [Xs Ys];
alldata = sortrows(alldata);
alldataxs = alldata(:,1);
alldatays = alldata(:,2);
% alldata = reshape(alldata, 100, []);
% alldata = mean(alldata, 1);
alldatays = movmean(alldatays, N/1e2);

figure();
hold on;
scatter(Xs, Ys);
x = -3:1e-2:-1;
plot(x, (3+x).^2/12, 'k');
x = -1:1e-2:1;
plot(x, 1/3*ones(size(x)), 'k');
x = 1:1e-2:3;
plot(x, (3-x).^2/12, 'k');
plot(alldataxs, alldatays, 'r');
% plot(linspace(-3, 3, length(alldata)), alldata, 'r');

MSE = mean((YYhat(:,1) - YYhat(:,2)).^2)

%% 1b
Ylin = X/5;

LMSE = mean((Ylin - Y).^2)

%% 2
muy = 1;
sigy = 8;
sigr = 5;
R = 10;

Ygaus = muy + sigy * randn(N, 1);
Rgaus = sigr * randn(N, R);

Xgaus = Ygaus + Rgaus;

% solve matrix equation Ax=b
% general case: 8.70; less general case: 8.77
% 8.77

% sig^2_X1 = var(y+r) = var(y) + var(r) = sigy^2 + sigr^2
% sig_X1X2 = E[(X1-mu_x1)*(X2-mu_x2)] = E[X1X2] - muY^2 = E[Y^2]-mu_Y^2 = var(Y)

% E[X1x2]=E[(Y+R1)(Y+R2)]=E[Y^2+YR_1+Y_R2+R_1R_2]=E[Y^2]

% diagonal terms are all : sigy^2+sigr^2
% all other terms are : sigy^2

% XiXj i != j

% A = [sigy^2+sigr^2 sigy^2; sigy^2 sigy^2+sigr^2]; 
% b = sigy^2 * ones(R, 1);
% a = A^-1 * b;

% 8.70
C_XX = sigy^2 * ones(R) + sigr^2 * eye(R);
C_XY = sigy^2 * ones(R, 1);
a = C_XX^-1 * C_XY;

% 8.74: comes from requirement that Yhat_lin is unbiased
% (same as 8.60 assuming that mu_x = mu_y since mu_r=0)
a0 = muy * (1 - sum(a));

Yhat = a0 + sum(a.' .* Xgaus, 2);
MSE = mean((Yhat - Ygaus).^2)

% 8.73: associated MMSE
sigy^2 - C_XY.'*(C_XX)^-1*C_XY

