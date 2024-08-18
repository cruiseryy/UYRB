clear;clc;close all 

addpath('utils_matlab');
addpath('data_matlab');

% read DPS policies using perfect inflow forecasts
tdata0 = importdata('re0.reference');
policy0 = tdata0(:,1:end-4);
obj0 = [];

for k = 1:size(tdata0,1)
    [hobj0, ~, ~, ~, ~] = his_func(policy0(k,:), 0);
    obj0(k,:) = mean(hobj0);
    k
end

% read DPS policies using binary inflow forecasts
tdata1 = importdata('re1.reference');
policy1 = tdata1(:,1:end-4);
obj1 = [];

for k = 1:size(tdata1,1)
    [hobj1, ~, ~, ~, ~] = his_func(policy1(k,:), 1);
    obj1(k,:) = mean(hobj1);
    k
end

% read the baseline performance using actual operations
load('his_obj.mat')
hobj = mean(his_obj);

% write the objs for later visualization
wobj0 = obj0 .* [-1, -1, 1, -1];
wobj1 = obj1 .* [-1, -1, 1, -1];
whobj = hobj .* [-1, -1, 1, -1];
writematrix(wobj0, 'figures/his_obj_dps0.csv')
writematrix(wobj1, 'figures/his_obj_dps1.csv')
writematrix(whobj, 'figures/his_obj_actual.csv')

% flip obj3 (water deficit) to align the preference directions across objs
obj0(:,3) = 1 - obj0(:,3);
obj1(:,3) = 1 - obj1(:,3); 
hobj(:,3) = 1 - hobj(:,3);

ymin = min([obj0; obj1; hobj]);
ymax = max([obj0; obj1; hobj]);

obj0 = (obj0 - repmat(ymin,size(obj0,1),1))./repmat(ymax-ymin,size(obj0,1),1);
obj1 = (obj1 - repmat(ymin,size(obj1,1),1))./repmat(ymax-ymin,size(obj1,1),1);
hobj = (hobj - ymin)./(ymax - ymin); 

% find the closest policy in the obj space measured by Euclidean distance
dev0 = sum((obj0 - repmat(hobj,size(obj0,1),1)).^2,2);
idx0 = find(dev0 == min(dev0));

dev1 = sum((obj1 - repmat(hobj,size(obj1,1),1)).^2,2);
idx1 = find(dev1 == min(dev1));

% a very rough parallel coordinate plot to check in the obj space
% flip all objs so that lower = better performance (a downward preference direction)
hobj = -hobj;
obj0 = -obj0;
obj1 = -obj1;
figure()
hold on
plot(1:4, hobj, 'kd-','LineWidth', 1.5, 'MarkerSize', 10);
plot(1:4, obj0(idx0,:) , 'rd-','LineWidth', 1.5, 'MarkerSize', 10)
plot(1:4, obj1(idx1,:) , 'bd-','LineWidth', 1.5, 'MarkerSize', 10)

for k = 1:size(tdata0,1)
    tplot = plot(1:4, obj0(k,:), 'r-','LineWidth',1);
    tplot.Color(4) = 0.15;
end

for k = 1:size(tdata1,1)
    tplot = plot(1:4, obj1(k,:), 'b-','LineWidth',1);
    tplot.Color(4) = 0.15;
end

plot([1 1]*1, [-0 1],'k')
plot([1 1]*2, [-0 1],'k')
plot([1 1]*3, [-0 1],'k')
plot([1 1]*4, [-0 1],'k')
axis([1 4 -1 0])
set(gca,'XTick',1:4)
set(gca,'XTickLabel', {'tot E (J1)','min E (J2)','WD (J3)','SP (J4)'})

legend('actual ops', 'DPS perfect', 'DPS binary')
grid

% get the operation trajs of DPS policies most similar to actual ops 
[idx0, idx1]
[~, dps_lyx_o0, dps_lyx_s0, dps_ljx_o0, dps_ljx_s0] = his_func(policy0(idx0,:), 0);
[~, dps_lyx_o1, dps_lyx_s1, dps_ljx_o1, dps_ljx_s1] = his_func(policy1(idx1,:), 1);

% get other commonly-used 'representative' policies
idx_others0 = [];
idx_others1 = [];
for k = 1:4
    idx_others0(k) = find(obj0(:,k) == min(obj0(:,k)));
    idx_others1(k) = find(obj1(:,k) == min(obj1(:,k)));
end

% to save the trajectories for plotting using python
lyx_outflow0 = [];
ljx_outflow0 = [];
lyx_storage0 = []; 
ljx_storage0 = [];
lyx_outflow1 = [];
ljx_outflow1 = [];
lyx_storage1 = []; 
ljx_storage1 = [];

% compare the trajs
load('flow_data_new.mat')
figure()
subplot(3,1,1)
hold on
cc_lyx_o = [];
for k = 1:4
    [~, tmp_o0, ~, ~, ~] = his_func(policy0(idx_others0(k),:), 0);
    [~, tmp_o1, ~, ~, ~] = his_func(policy1(idx_others1(k),:), 1);
    tplot = plot(tmp_o0, 'r-','LineWidth',2);
    tplot.Color(4) = 0.15;
    tplot = plot(tmp_o1, 'b-','LineWidth',2);
    tplot.Color(4) = 0.15;
    [k, corr(tmp_o0,lyx_out), corr(tmp_o1,lyx_out)]
    
    lyx_outflow0 = [lyx_outflow0, tmp_o0];
    lyx_outflow1 = [lyx_outflow1, tmp_o1];
end 
plot(dps_lyx_o0, 'r-','LineWidth',2)
plot(dps_lyx_o1, 'b-','LineWidth',2)
plot(lyx_out, 'k-', 'LineWidth',2)
[5, corr(dps_lyx_o0,lyx_out), corr(dps_lyx_o1,lyx_out)]
lyx_outflow0 = [lyx_outflow0, dps_lyx_o0, lyx_out];
lyx_outflow1 = [lyx_outflow1, dps_lyx_o1, lyx_out];

subplot(3,1,2)
hold on
cc_ljx_o = [];
for k = 1:4
    [~, ~, ~, tmp_o0, ~] = his_func(policy0(idx_others0(k),:), 0);
    [~, ~, ~, tmp_o1, ~] = his_func(policy1(idx_others1(k),:), 1);
    tplot = plot(tmp_o0, 'r-','LineWidth',2);
    tplot.Color(4) = 0.15;
    tplot = plot(tmp_o1, 'b-','LineWidth',2);
    tplot.Color(4) = 0.15;
    [k, corr(tmp_o0,ljx_out), corr(tmp_o1,ljx_out)]

    ljx_outflow0 = [ljx_outflow0, tmp_o0];
    ljx_outflow1 = [ljx_outflow1, tmp_o1];
end 
plot(dps_ljx_o0, 'r-','LineWidth',2)
plot(dps_ljx_o1, 'b-','LineWidth',2)
plot(ljx_out, 'k-', 'LineWidth',2)
[5, corr(dps_ljx_o0,ljx_out), corr(dps_ljx_o1,ljx_out)]
ljx_outflow0 = [ljx_outflow0, dps_ljx_o0, ljx_out];
ljx_outflow1 = [ljx_outflow1, dps_ljx_o1, ljx_out];

subplot(3,1,3)
hold on
cc_s = [];
for k = 1:4
    [~, ~, tmp_s10, ~, tmp_s20] = his_func(policy0(idx_others0(k),:), 0);
    [~, ~, tmp_s11, ~, tmp_s21] = his_func(policy1(idx_others1(k),:), 1);
    tmp_s0 = tmp_s10 + tmp_s20; 
    tmp_s1 = tmp_s11 + tmp_s21;
    tplot = plot(tmp_s0, 'r-','LineWidth',2);
    tplot.Color(4) = 0.15;
    tplot = plot(tmp_s1, 'b-','LineWidth',2);
    tplot.Color(4) = 0.15;
    [k, corr(tmp_s0,lyx_s'+ljx_s'), corr(tmp_s1,lyx_s'+ljx_s')]

    lyx_storage0 = [lyx_storage0, tmp_s10];
    ljx_storage0 = [ljx_storage0, tmp_s20];
    lyx_storage1 = [lyx_storage1, tmp_s11];
    ljx_storage1 = [ljx_storage1, tmp_s21];
end 
plot(dps_lyx_s0 + dps_ljx_s0, 'r-','LineWidth',2)
plot(dps_lyx_s1 + dps_ljx_s1, 'b-','LineWidth',2)
plot(lyx_s + ljx_s, 'k-', 'LineWidth',2)
[5, corr(dps_lyx_s0+dps_ljx_s0,lyx_s'+ljx_s'), corr(dps_lyx_s1+dps_ljx_s1,lyx_s'+ljx_s')]
lyx_storage0 = [lyx_storage0, dps_lyx_s0, lyx_s'];
ljx_storage0 = [ljx_storage0, dps_ljx_s0, ljx_s'];
lyx_storage1 = [lyx_storage1, dps_lyx_s1, lyx_s'];
ljx_storage1 = [ljx_storage1, dps_ljx_s1, ljx_s'];

writematrix(lyx_outflow0, 'figures/lyx_outflow0.csv')
writematrix(lyx_outflow1, 'figures/lyx_outflow1.csv')
writematrix(ljx_outflow0, 'figures/ljx_outflow0.csv')
writematrix(ljx_outflow1, 'figures/ljx_outflow1.csv')
writematrix(lyx_storage0, 'figures/lyx_storage0.csv')
writematrix(lyx_storage1, 'figures/lyx_storage1.csv')
writematrix(ljx_storage0, 'figures/ljx_storage0.csv')
writematrix(ljx_storage1, 'figures/ljx_storage1.csv')


% The representative policies minimzing the 2nd obj share relatively 
% better similarity (measured by Pearson's CC) with the historical
% actual operations compared to other representative policies. But overall
% they are beaten by the preference-preserving policies selected with our
% framework.





