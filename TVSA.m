clear;clc;

addpath('utils_matlab');
addpath('data_matlab');

rng(0)

% load parameters
load('paras.mat');

% read DPS policies & obj values
tdata0 = importdata('re0.reference');
obj0 = tdata0(:,end-3:end);
policy0 = tdata0(:,1:end-4);
tdata1 = importdata('re1.reference');
policy1 = tdata1(:,1:end-4);
obj1 = tdata1(:,end-3:end);

% Representative policies are selected to mimic actual operations
% his_eval_dps.m returns the indices of selected rep policies
idx0 = 153;
idx1 = 56;

% read synthetic flow data
inf_pool = importdata('inflow.txt');
lat_pool = importdata('latflow.txt');

na = 5; nb = 5;
s1 = zeros(2,1); s2 = zeros(2,1);
Cp = zeros(1,na*nb);
Rp = zeros(1,na*nb);
Wp1 = zeros(1,na); Wp2 = zeros(1,na);

% select a DPS policy for 
type = 1;
if type == 0
    policy = policy0(idx0, :);
else
    policy = policy1(idx1, :);
end

% Parse policy into RBF coefficients
mark = 1;
for i = 1:na
    for j = 1:nb
        Cp(na*(i-1)+j) = policy(mark);
        mark = mark + 1;
    end
    for j = 1:nb
        Rp(na*(i-1)+j) = policy(mark);
        mark = mark + 1;
    end
    Wp1(i) = policy(mark);
    mark = mark + 1;
    Wp2(i) = policy(mark);
    mark = mark + 1;
end

% Normalize weights and scale radius 
Wp1 = Wp1/sum(Wp1);
Wp2 = Wp2/sum(Wp2);
Rp = 10*Rp;

% 100 samples (SOWs) are used to estimate the variances and their components
xx = zeros(100,5,36); % sampled-based variance estimates of inputs
yy = zeros(100,2,36); % sampled-based variance estimates of inputs
tmpobj = zeros(100,4);

for k = 1:100
    tmp_idx = randi(20000);
    tinf = inf_pool(tmp_idx,:);
    tlat = lat_pool(tmp_idx,:);
    
    s1(:) = 0;s1(1) = 120e8;
    s2(:) = 0;s2(1) = 25e8;
    
    tote = 0;
    mine = 1e8;
    water_deficit = 0;
    tmpmr = 0;
    
    for i = 1:36
        
        tmpinf = exp(tinf(i)*std_mm(i) + inf_mm(i));
        tmplat = b1 + b2*tmpinf + tlat(i)/4*(db1 + db2*tmpinf);
        tmpinf = tmpinf*10*24*3600;
        tmplat = tmplat*10*24*3600;
        
        if type == 0
            xvar = [(s1(i)-smin1)/(smax1-smin1) (s2(i)-smin2)/(smax2-smin2) mod(i-1,36)/36 (tinf(i)+2)/4 (tlat(i)+2)/4]; % ref0 xvar
        else
            xvar = [(s1(i)-smin1)/(smax1-smin1) (s2(i)-smin2)/(smax2-smin2) mod(i-1,36)/36 (tinf(i)>0)*0.5+0.25 (tlat(i)>0)*0.5+0.25]; % ref0 xvar
        end
        xx(k,:,i) = xvar;
        yy(k,:,i) = [RBF(xvar,Cp,Rp,Wp1) RBF(xvar,Cp,Rp,Wp2)]; 
        
        tmpd1 = RBF(xvar,Cp,Rp,Wp1)*(rmax1(i)-rmin1(i)) + rmin1(i);
        tmpd1 = tmpd1*10*24*3600;
        tmpa1 = max([0 s1(i)+tmpinf-smax1-tmpd1]); 
        tmpr1 = min([tmpa1+tmpd1 s1(i)+tmpinf-smin1]);
        s1(i+1) = s1(i) + tmpinf - tmpr1;
        
        tmpd2 = RBF(xvar,Cp,Rp,Wp2)*(rmax2(i)-rmin2(i)) + rmin2(i);
        tmpd2 = tmpd2*10*24*3600;
        tmpa2 = max([0 s2(i)+tmplat+tmpr1-smax2-tmpd2]);
        tmpr2 = min([tmpd2+tmpa2 s2(i)+tmplat+tmpr1-smin2]);
        s2(i+1) = s2(i) + tmpr1 + tmplat - tmpr2;
        
        tmp_hp = 0;
        tmp_hp = tmp_hp + 8.5*1e3*min([tmpr1 tmpd1])*(inhead1(s1(i))-outhead1(tmpr1));
        tmp_hp = tmp_hp + 8.5*1e3*min([tmpr2 tmpd2])*(inhead2(s2(i))-outhead2(tmpr2));
        for kk = 1:10
            tmp_hp = tmp_hp + 8.5*1e3*kkhead(kk)*(tmpr1 + tmplat*kkcsa(kk));
        end
        tmp_hp = tmp_hp/1e8/1e3/3600;
        
        mine = min(tmp_hp, mine);
        tote = tote + tmp_hp; 
        
        tmpmr = tmpmr + tmpr2/10/24/3600/3;
        if mod(i,3) == 0 
            midx = floor((i-0.1)/3) + 1;
            water_deficit = water_deficit + 1/36 * max(0, 1.2*demand(midx) - tmpmr) / 1.2 / demand(midx);
            tmpmr = 0;
        end     
    end
    storage_penalty = (s1(end) + s2(end))/(smax1 + smax2);
    tmpobj(k,:) = [tote mine water_deficit storage_penalty];
    
    
end

% for both row and column indexing:
% 1: lyx storage; 2: ljx storage; 3: T-o-Y; 4: inflow; 5: lateral flow
% for example, (1, 2, 5) is the var comp of interaction between lyx storage
% and ljx storage in May's decisions
var_comp_lyx = zeros(5, 5, 36);
var_comp_ljx = zeros(5, 5, 36);

for k = 1:36
    for m = 1:5
        for l = 1:5
            tmp_cov = cov(xx(:,m,k),xx(:,l,k));
            tmp_cov = tmp_cov(2);
            for si = 1:100
                var_comp_lyx(m,l,k) = var_comp_lyx(m,l,k) + 1/100*RBF_derivative(xx(i,:,k),Cp,Rp,Wp1,m)*RBF_derivative(xx(i,:,k),Cp,Rp,Wp1,l)*tmp_cov;
                var_comp_ljx(m,l,k) = var_comp_ljx(m,l,k) + 1/100*RBF_derivative(xx(i,:,k),Cp,Rp,Wp2,m)*RBF_derivative(xx(i,:,k),Cp,Rp,Wp2,l)*tmp_cov;
            end
        end
    end
end

% uncomment this chunk of code to evaluate the TVSA decomposition
% sample-based total variances of Y are compared againts the estimates from
% summing up all contributing components from Xs.
% ONLY NOT statistically significant for LJX release decisions with perfect
% inflow data. This is largely due to ONE outlier (index = 20), and R^2 is
% improved to > 0.8 when the outlier is removed. 
yy1 = reshape(yy(:,1,:),100,36);
yyvar1 = squeeze(var(yy1,[],1));
re_yyvar1 = squeeze(sum(var_comp_lyx, [1,2]));
[b1,~,~,~,stats1] = regress(yyvar1',[ones(36,1), re_yyvar1])
yy2 = reshape(yy(:,2,:),100,36);
yyvar2 = squeeze(var(yy2,[],1));
re_yyvar2 = squeeze(sum(var_comp_ljx, [1,2]));
[b2,~,~,~,stats2] = regress(yyvar2',[ones(36,1), re_yyvar2])
figure(1)
subplot(1,2,1)
scatter(yyvar1, re_yyvar1, 'rd')
subplot(1,2,2)
scatter(yyvar2, re_yyvar2, 'rd')
if type == 0
    hold on
    scatter(yyvar2(20), re_yyvar2(20), 'black','filled','o')
    yyvar2(20) = [];
    re_yyvar2(20) = [];
    [b2,~,~,~,stats2] = regress(yyvar2',[ones(35,1), re_yyvar2])
end

pause = 1;

lyx_var = zeros(5, 36);
ljx_var = zeros(5, 36);

tmark = 1;
for k = 1:4
    lyx_var(tmark,:) = var_comp_lyx(tmark,tmark,:);
    ljx_var(tmark,:) = var_comp_ljx(tmark,tmark,:);
    tmark = tmark + 1;
    % skip the ToY
    if tmark == 3
        tmark = tmark + 1;
    end
end

lyx_var_tot = squeeze(sum(var_comp_lyx, [1,2]));
ljx_var_tot = squeeze(sum(var_comp_ljx, [1,2]));

lyx_var(5,:) = lyx_var_tot' - sum(lyx_var(1:4,:));
ljx_var(5,:) = ljx_var_tot' - sum(ljx_var(1:4,:));

lyx_inter = zeros(5, 36);
ljx_inter = zeros(5, 36);

lyx_inter(1,:) = var_comp_lyx(4,5,:); % Inf - Lat
lyx_inter(2,:) = var_comp_lyx(4,1,:); % Inf - LYX storage
lyx_inter(3,:) = var_comp_lyx(4,2,:); % Inf - LJX storage
lyx_inter(4,:) = var_comp_lyx(5,1,:); % Lat - LYX storage
lyx_inter(5,:) = var_comp_lyx(5,2,:); % Lat - LJX storage

ljx_inter(1,:) = var_comp_ljx(4,5,:); % Inf - Lat
ljx_inter(2,:) = var_comp_ljx(4,1,:); % Inf - LYX storage
ljx_inter(3,:) = var_comp_ljx(4,2,:); % Inf - LJX storage
ljx_inter(4,:) = var_comp_ljx(5,1,:); % Lat - LYX storage
ljx_inter(5,:) = var_comp_ljx(5,2,:); % Lat - LJX storage

% a simple visualization to examine the decomposition results
tmpvar = lyx_inter;
figure(11)
subplot(2,1,type+1)
hold on 
for k = 1:36
    tmpmark = 0;
    
    tplot = fill(k+[-0.25 -0.25 0.25 0.25 -0.25],tmpmark+[0 tmpvar(1,k) tmpvar(1,k) 0 0],'g','EdgeColor','g'); % S1*S1
    tplor.Color(4) = 0.4;
    tmpmark = tmpmark + tmpvar(1,k);
    
    tplot = fill(k+[-0.25 -0.25 0.25 0.25 -0.25],tmpmark+[0 tmpvar(2,k) tmpvar(2,k) 0 0],'m','EdgeColor','m'); % S2*S2
    tplor.Color(4) = 0.4;
    tmpmark = tmpmark + tmpvar(2,k);
    
    tplot = fill(k+[-0.25 -0.25 0.25 0.25 -0.25],tmpmark+[0 tmpvar(3,k) tmpvar(3,k) 0 0],'r','EdgeColor','r'); % I*I
    tplor.Color(4) = 0.4;
    tmpmark = tmpmark + tmpvar(3,k);
    
    tplot = fill(k+[-0.25 -0.25 0.25 0.25 -0.25],tmpmark+[0 tmpvar(4,k) tmpvar(4,k) 0 0],'b','EdgeColor','b'); % L*L
    tplor.Color(4) = 0.4;
    tmpmark = tmpmark + tmpvar(4,k);
    
    plot(k+[-0.25 -0.25 0.25 0.25 -0.25],tmpmark + [0 tmpvar(5,k) tmpvar(5,k) 0 0],'k--')  %interaction
    tmpmark = tmpmark + tmpvar(5,k);

    [k, tmpvar(2,k)/tmpmark]
end 
grid


tmpvar = ljx_inter;
figure(12)
subplot(2,1,type+1)
hold on 
for k = 1:36
    tmpmark = 0;
    
    tplot = fill(k+[-0.25 -0.25 0.25 0.25 -0.25],tmpmark+[0 tmpvar(1,k) tmpvar(1,k) 0 0],'g','EdgeColor','g');
    tplor.Color(4) = 0.4;
    tmpmark = tmpmark + tmpvar(1,k);
    
    tplot = fill(k+[-0.25 -0.25 0.25 0.25 -0.25],tmpmark+[0 tmpvar(2,k) tmpvar(2,k) 0 0],'m','EdgeColor','m');
    tplor.Color(4) = 0.4;
    tmpmark = tmpmark + tmpvar(2,k);
    
    tplot = fill(k+[-0.25 -0.25 0.25 0.25 -0.25],tmpmark+[0 tmpvar(3,k) tmpvar(3,k) 0 0],'r','EdgeColor','r');
    tplor.Color(4) = 0.4;
    tmpmark = tmpmark + tmpvar(3,k);
    
    tplot = fill(k+[-0.25 -0.25 0.25 0.25 -0.25],tmpmark+[0 tmpvar(4,k) tmpvar(4,k) 0 0],'b','EdgeColor','b');
    tplor.Color(4) = 0.4;
    tmpmark = tmpmark + tmpvar(4,k);
    
    plot(k+[-0.25 -0.25 0.25 0.25 -0.25],tmpmark + [0 tmpvar(5,k) tmpvar(5,k) 0 0],'k--')  %interaction
    tmpmark = tmpmark + tmpvar(5,k);

    [k, tmpvar(2,k)/tmpmark]
end 
grid














