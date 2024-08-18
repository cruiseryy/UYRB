clear;clc;close all

addpath('utils_matlab');
addpath('data_matlab');

rng(0)

% read DPS policies & obj values
tdata0 = importdata('re0.reference');
obj0 = tdata0(:,end-3:end);
policy0 = tdata0(:,1:end-4);
tdata1 = importdata('re1.reference');
policy1 = tdata1(:,1:end-4);
obj1 = tdata1(:,end-3:end);

% create a set of altered flow conditions
% the amplification factors range from [1, 5] and are sampled using LHS
n = 1000;
af = 1 + repmat([4 4],n,1) .* lhsdesign(n,2); 

% Representative policies are selected to mimic actual operations
% his_eval_dps.m returns the indices of selected rep policies
idx0 = 153;
idx1 = 56;

% MARK
% base objective
obj_all = [obj0; obj1];
obj_base = quantile(obj_all, 0.9) .* [-1, -1, 1, -1];

reliability = zeros(n,2); 

for k = 1:n
    tic 
    tmpr0 = robust_func(policy0(idx0,:), 0, af(k,:), obj_base);
    tmpr1 = robust_func(policy1(idx1,:), 1, af(k,:), obj_base);
    reliability(k,:) = [tmpr0 tmpr1];
    [k toc]
end

% simple visualization to check results
cmin = min(reliability(:));
cmax = max(reliability(:));
base = [1 1];
figure()
subplot(1,3,1)
hold on
scatter(af(:,1), af(:,2), [], reliability(:,1), 'filled')
scatter(base(1), base(2), 100, 'kx')
th = colorbar;
caxis([cmin, cmax]);
ylabel(th, 'Reliability (%)')
grid
xlabel('inflow amplification factor')
ylabel('latflow amplification factor')
title('(a) True')
subplot(1,3,2)
hold on
scatter(af(:,1), af(:,2), [], reliability(:,2), 'filled')
scatter(base(1), base(2), 100, 'kx')
th = colorbar;
caxis([cmin, cmax]);
ylabel(th, 'Reliability (%)')
grid
xlabel('inflow amplification factor')
ylabel('latflow amplification factor')
title('(b) Binary')
subplot(1,3,3)
hold on
scatter(af(:,1), af(:,2), [], reliability(:,1) - reliability(:,2), 'filled')
scatter(base(1), base(2), 100, 'kx')
th = colorbar;
ylabel(th, 'Reliability (%)')
grid
xlabel('inflow amplification factor')
ylabel('latflow amplification factor')
title('(c) Difference')

[corr(reliability(:,1),af(:,1)), corr(reliability(:,1),af(:,2))]
[corr(reliability(:,2),af(:,1)), corr(reliability(:,2),af(:,2))]
[corr(reliability(:,1)-reliability(:,2),af(:,1)), corr(reliability(:,1)-reliability(:,2),af(:,2))]

writematrix(af, 'figures/LHS_af.csv')
writematrix(reliability, 'figures/robustness.csv')

% the dps1 reliability as well as the reliability diff are negatively and
% positively correlated with latflow af, both being satistically
% significant while the dps0 reliability show no correlation with the
% latflow af. Note that latflow only plays a very small portion here!!!
% also, the reliability diff highly positively correlated with inflow AF,
% indicating that in more altered worlds, the high-quality inflow forecasts
% play more important roles!
