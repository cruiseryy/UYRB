function [obj, r1, s1, r2, s2] = his_func(policy, ref_type)

% load parameters
load('paras.mat');

% read historic 10-day flow data 2001-2009 
load('flow_data_new.mat');


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

s1 = zeros(2,1);
s2 = zeros(2,1);

% initial reservoir storage
s1(1) = 120e8;
% s1(1) = 112e8; % pre-dam LYX water head = 2557m in 2001 Jan
s2(1) = 25e8;
% s2(1) = 22e8; % pre-dam LJX water head = 1719 m in 2001 Jan
r1 = zeros(2,1); r2 = zeros(2,1);

tmark = 1;
hydro_power = [];
storage_penalty = [];
tmpmr = 0; 
water_deficit = zeros(9,1);

% Simulation loop for 9 years
for y = 1:9
    for i = 1:36
        
        tmpinf = lyx_in(tmark); tmplat = ljx_in(tmark) - lyx_out(tmark); 
        infx = (log(tmpinf) - inf_mm(i))/std_mm(i);
        lafx = (tmplat - b1 - b2*tmpinf)*4/(db1 + db2*tmpinf);
        tmpinf = tmpinf*10*24*3600;
        tmplat = tmplat*10*24*3600;
        
        if ref_type == 0
            xvar = [(s1(tmark)-smin1)/(smax1-smin1) (s2(tmark)-smin2)/(smax2-smin2) mod(i-1,36)/36 (infx+2)/4 (lafx+2)/4]; % ref0 xvar
        else
            xvar = [(s1(tmark)-smin1)/(smax1-smin1) (s2(tmark)-smin2)/(smax2-smin2) mod(i-1,36)/36 (infx>0)*0.5+0.25 (lafx>0)*0.5+0.25]; % ref1 xvar
        end

        tmpd1 = RBF(xvar,Cp,Rp,Wp1) * (rmax1(i)-rmin1(i)) + rmin1(i);
        tmpd1 = tmpd1*10*24*3600;
        tmpa1 = max([0 s1(tmark)+tmpinf-smax1-tmpd1]); % estimate spillage to ensure no overflow
        tmpr1 = min([tmpa1+tmpd1 s1(tmark)+tmpinf-smin1]); % estimate actual release and ensure no min storage level
        s1(tmark+1) = s1(tmark) + tmpinf - tmpr1;
        r1(tmark) = tmpr1/10/24/3600;
        
        tmpd2 = RBF(xvar,Cp,Rp,Wp2) * (rmax2(i)-rmin2(i)) + rmin2(i);
        tmpd2 = tmpd2*10*24*3600;
        tmpa2 = max([0 s2(tmark)+tmplat+tmpr1-smax2-tmpd2]);
        tmpr2 = min([tmpd2+tmpa2 s2(tmark)+tmplat+tmpr1-smin2]);
        s2(tmark+1) = s2(tmark) + tmpr1 + tmplat - tmpr2;
        r2(tmark) = tmpr2/10/24/3600;
        
        tmp_hp = 0; 
        tmp_hp = tmp_hp + 8.5*1e3*min([tmpr1 tmpd1])*(inhead1(s1(tmark))-outhead1(tmpr1));
        tmp_hp = tmp_hp + 8.5*1e3*min([tmpr2 tmpd2])*(inhead2(s2(tmark))-outhead2(tmpr2));
        for k = 1:10
            tmp_hp = tmp_hp + 8.5*1e3*kkhead(k)*(tmpr1 + tmplat*kkcsa(k));
        end
        tmp_hp = tmp_hp/1e8/1e3/3600;
        hydro_power(tmark) = tmp_hp;
        
        tmpmr = tmpmr + tmpr2/10/24/3600/3;
        if mod(i,3) == 0
            midx = floor((i-0.1)/3) + 1;
            water_deficit(y) = water_deficit(y) + 1/36 * max(0, 1.2*demand(midx) - tmpmr) / 1.2 / demand(midx);
            tmpmr = 0;
        end
        tmark = tmark + 1;
    end
    storage_penalty(y) = (s1(tmark) + s2(tmark))/(smax1+smax2);
end

tote = [];
mine = [];
for y = 1:9
    tmp_hp = hydro_power((y-1)*36+1:36*y);
    tote(y) = sum(tmp_hp);
    mine(y) = min(tmp_hp);
end

obj = [tote' mine' water_deficit storage_penalty'];

end




