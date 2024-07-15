clear;clc;close all
rng(2)
smin1 = 50e8;
smax1 = 247e8;
smin2 = 6e8;
smax2 = 57e8;

tdata0 = importdata('../../re0.reference');
policy0 = tdata0(:,1:end-4);
tdata1 = importdata('../../re1.reference');
policy1 = tdata1(:,1:end-4);

n = 100;
ss = lhsdesign(n,2);
s0 = 1 + repmat([4 4],n,1).*ss;

ref0_idx = 165;
ref1_idx = 56;

% obj_base = [324.4417    5.8861    0.0727    0.5763]; % sum(diff>0) = 96
obj_base = [220 4 0.11 0.4]; % sum(diff>0) = 100

reliability = zeros(n,2); 

for k = 1:n
    tic 
    tmpr0 = rb_flow_new(policy0(ref0_idx,:), 0, s0(k,:), obj_base);
    tmpr1 = rb_flow_new(policy1(ref1_idx,:), 1, s0(k,:), obj_base);
    reliability(k,:) = [tmpr0 tmpr1];
    [k toc]
   
end
% ts = [5*ones(50,1) linspace(1,5,50)'];
% tr0 = zeros(50,1); tr1 = zeros(50,1);
% for k = 1:50
%     tr0(k) = rb_flow_new(policy0(ref0_idx,:), 0, ts(k,:), obj_base);
%     tr1(k) = rb_flow_new(policy1(ref1_idx,:), 1, ts(k,:), obj_base);
% end
% figure()
% subplot(1,2,1)
% scatter(ts(:,2), tr0)
% subplot(1,2,2)
% scatter(ts(:,2), tr1)

sbase = [1 1];
figure()
subplot(1,3,1)
hold on
scatter(s0(:,1), s0(:,2), [], reliability(:,1), 'filled')
scatter(sbase(1), sbase(2), 100, 'kx')
th = colorbar;
ylabel(th, 'Reliability (%)')
axis([1 5 1 5])
grid
xlabel('inflow rescale coef')
ylabel('latflow rescale coef')
title('(a) True')
subplot(1,3,2)
hold on
scatter(s0(:,1), s0(:,2), [], reliability(:,2), 'filled')
scatter(sbase(1), sbase(2), 100, 'kx')
th = colorbar;
ylabel(th, 'Reliability (%)')
axis([1 5 1 5])
grid
xlabel('inflow rescale coef')
ylabel('latflow rescale coef')
title('(b) Binary')
subplot(1,3,3)
hold on
scatter(s0(:,1), s0(:,2), [], reliability(:,1) - reliability(:,2), 'filled')
scatter(sbase(1), sbase(2), 100, 'kx')
th = colorbar;
ylabel(th, 'Reliability (%)')
grid
axis([1 5 1 5])
xlabel('inflow rescale coef')
ylabel('latflow rescale coef')
title('(c) Difference')
