% clc; clear; close all;

%% 
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

