%% Modeling: Clustered network
clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\modeling_3\
addpath(genpath(pwd))
addpath C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw\
addpath C:\Users\duodenum\Desktop\brain_stuff\misc\2019_03_03_BCT
%% Parameters
b = -3; 
tau = 0.1;
C = 1;
s = 0; % rest
s_task = 3;
k = 1;
tspan = 300;
dt = 10 / 1000;

W = load('C:\Users\duodenum\Desktop\brain_stuff\modelingresources\averageConnectivity_Fpt.mat');
W = 10.^(W.Fpt);
n = length(W);
W(logical(eye(n))) = 1;
W(~logical(eye(n))) = 2*(W(~logical(eye(n))) ./ sum(W(~logical(eye(n))), 'all'));

nsims = 30;
ntp = 20;
acw0s_rest = zeros(n, ntp, nsims);
acw0s_task = zeros(n, ntp, nsims);
means_rest = zeros(n, nsims);
means_task = zeros(n, nsims);
vars_rest = zeros(n, nsims);
vars_task = zeros(n, nsims);

parfor i = 1:nsims
    tic
    %% Set up W
    %% run simulations
    x = wilsoncowan_RK2(tau, b, W, k, s, C, tspan);
    acw0s_rest(:, :, i) = acw_windowed_loop(x, 1/dt, 10);
    means_rest(:, i) = mean(x, 2);
    var_rest(:, i) = std(x, [], 2);

    x = wilsoncowan_RK2(tau, b, W, k, s_task, C, tspan);
    acw0s_task(:, :, i) = acw_windowed_loop(x, 1/dt, 10);
    means_task(:, i) = mean(x, 2);
    var_task(:, i) = std(x, [], 2);
    toc
    fprintf("\ni = %d\n", i)
end

save("dtinetwork.mat", "acw0s_rest", "acw0s_task", "means_rest", "means_task", ...
    "var_task", "var_rest")
