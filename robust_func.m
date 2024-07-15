function count = robust_func(policy, type, alpha, obj_base)

% load parameters
load('paras.mat');

% read synthetic flow data
inf_pool = importdata('inflow.txt');

na = 5; nb = 5;
Cp = zeros(1,na*nb);
Rp = zeros(1,na*nb);
Wp1 = zeros(1,na); Wp2 = zeros(1,na);

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

s1 = zeros(37,1);
s2 = zeros(37,1);
tmpobj = zeros(2,4);

% generate altered SOWs based on the scaling coef by re-sampling from the og synthetic inflow data
SOW = sow_idx(alpha(1),100);
latidx = alter_normal(alpha(2),100,36);

for k = 1:100
    % initialize the inflow data & storage state
    tinf = inf_pool(SOW(k),:);
    tlat = latidx(k,:);
    s1(:) = 0; s1(1) = 120e8; 
    s2(:) = 0; s2(1) = 25e8;
    
    % reset the obj values
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

count = 0;
for k = 1:100
    count = count + (tmpobj(k,1)>obj_base(1))*(tmpobj(k,2)>obj_base(2))*(tmpobj(k,3)<obj_base(3))*(tmpobj(k,4)>obj_base(4));
end

end














