clear;clc;close all
rng(0)

load('flow_data_new.mat')

y0 = reshape(lyx_in,36,9)';

y1 = log(y0);
mm_mean = mean(y1);
mm_std = std(y1,[],1);

ya = (y1 - repmat(mm_mean,9,1))./repmat(mm_std,9,1);

% after step 1 log() to remove skewness
% after step 2 (y - mean(y))/std(y) for standardization (to a roughly
% normal distribution)
% figure()
% boxplot(ya)
% xlabel('Time')
% ylabel('Y anomaly')
% set(gca,'XTick',2:3:36)
% set(gca,'XTickLabel',1:12)
% grid

% n = 9e2; 
n = 2e4;
% m = randi(9,n,36);
% x = ya(m);
x = randn(n,36);

ccy = corr(ya);
q = chol(ccy + 1e-4*diag(ones(36,1)));

z = x*q;

latz = randn(n, 36);

dlmwrite('inflow.txt',z,'delimiter',' ')
dlmwrite('latflow.txt',latz,'delimiter',' ')

% ccz = corr(z); 

% figure()
% subplot(1,2,1)
% pcolor(ccy)
% h = colorbar;
% ylabel(h,'CC')
% xlabel('lag')
% ylabel('lag')
% title('corr(y)')
% 
% subplot(1,2,2)
% pcolor(ccz)
% h = colorbar;
% ylabel(h,'CC')
% xlabel('lag')
% ylabel('lag')
% title('corr(z)')
% 
% figure()
% 
% subplot(2,1,1)
% hold on 
% 
% cc1lag_y = zeros(1,35);
% cc1lag_z = zeros(100,35);
% 
% for k = 1:35
%     cc1lag_y(k) = corr(ya(:,k), ya(:,k+1)); 
%     for i = 1:100
%         cc1lag_z(i,k) = corr(z(9*(i-1)+1:9*i, k), z(9*(i-1)+1:9*i, k+1));
%     end
%     plot(k*[1 1], [quantile(cc1lag_z(:,k),0.1) quantile(cc1lag_z(:,k),0.9)], 'b-')
%     plot(k, cc1lag_y(k), 'rd')
% end
% title('time-varying lag-1 CC')
% ylabel('CC')
% xlabel('T_0')
% grid
% 
% subplot(2,1,2)
% hold on 
% 
% cc2lag_y = zeros(1,34);
% cc2lag_z = zeros(100,34);
% 
% for k = 1:34
%     cc1lag_y(k) = corr(ya(:,k), ya(:,k+2)); 
%     for i = 1:100
%         cc1lag_z(i,k) = corr(z(9*(i-1)+1:9*i, k), z(9*(i-1)+1:9*i, k+2));
%     end
%     plot(k*[1 1], [quantile(cc1lag_z(:,k),0.1) quantile(cc1lag_z(:,k),0.9)], 'b-')
%     plot(k, cc1lag_y(k), 'rd')
% end
% title('time-varying lag-2 CC')
% ylabel('CC')
% xlabel('T_0')
% grid












