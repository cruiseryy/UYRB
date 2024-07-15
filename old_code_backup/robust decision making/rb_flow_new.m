function count = rb_flow_new(policy, type, alpha, obj_base)
% rng(0)

% read synthetically generated flow data 
inf_pool = importdata('../streamflow_generator/inflow.txt');
lat_pool = importdata('../streamflow_generator/latflow.txt');

% read some parameters
b1 = 12.8416;
b2 = 0.2515;
db1 = 35.824;
db2 = 0.0521;
smin1 = 50e8;
smax1 = 247e8;
smin2 = 6e8;
smax2 = 57e8;
kkhead = [205,14,122,12.5,18.7,99.3,16,16,73,9.2];
kkcsa = [0.0119,0.0151,0.1062,0.1112,0.1310,0.2426,0.2651,0.2849,0.3047,0.3262];
rmin1 = [285.6,238.4,169.6,142.4,133.6,149.6,167.2,168.0,162.4,196.8,332.8,356.8,293.6,208.0,366.4,299.2,263.2,237.6,316.0,305.6,363.2,239.2,189.6,191.2,203.2,277.6,275.2,291.2,228.8,249.6,325.6,320.8,276.0,291.2,274.4,269.6];
rmax1 = [745.2,615.6,588.0,628.8,662.4,750.0,668.4,776.4,898.8,1046.4,1136.4,1162.8,1123.2,1160.4,1059.6,1053.6,1050.0,1039.2,1009.2,1034.4,993.6,847.2,926.4,886.8,860.4,747.6,866.4,700.8,922.8,962.4,914.4,871.2,783.6,753.6,738.0,672.0];
rmin2 = [238.4,233.6,231.2,223.2,208.0,191.2,186.4,182.4,189.6,226.4,234.4,379.2,729.6,405.6,277.6,500.8,450.4,485.6,461.6,376.0,295.2,298.4,315.2,374.4,364.0,330.4,413.6,429.6,536.8,576.8,615.2,457.6,336.0,345.6,328.0,292.0];
rmax2 = [555.6,552.0,554.4,546.0,537.6,445.2,372.0,718.8,1326.0,1282.8,1380.0,1358.4,1413.6,1377.6,1286.4,1375.2,1407.6,1375.2,1260.0,1140.0,1195.2,1099.2,1099.2,1243.2,1036.8,1095.6,1142.4,1358.4,1423.2,1549.2,1558.8,960.0,678.0,596.4,601.2,602.4];
demand = [550,500,350,750,1100,900,800,750,750,800,650,600];
inf_mm = [5.07559,4.98739,5.02828,5.03286,5.09163,5.17904,5.22714,5.36965,5.43714,5.57579,5.64422,5.87030,5.90242,5.94810,6.14867,6.38404,6.54943,6.78031,6.87444,6.90823,6.80754,6.72528,6.76357,6.78094,6.86762,6.78863,6.83915,6.88035,6.72561,6.49341,6.23301,6.01686,5.73047,5.51399,5.35356,5.21524];
std_mm = [0.34692,0.30352,0.23062,0.20955,0.22309,0.20396,0.16340,0.16398,0.19826,0.16088,0.21378,0.30881,0.37468,0.32197,0.26952,0.27893,0.32660,0.39741,0.31701,0.45362,0.45336,0.42666,0.41728,0.47373,0.50394,0.41208,0.35528,0.46297,0.45632,0.38720,0.32151,0.29521,0.28154,0.24770,0.26209,0.22290];
na = 5; nb = 5;
Cp = zeros(1,na*nb);
Rp = zeros(1,na*nb);
Wp1 = zeros(1,na); Wp2 = zeros(1,na);

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

Wp1 = Wp1/sum(Wp1);
Wp2 = Wp2/sum(Wp2);
Rp = 10*Rp;

s1 = zeros(37,1);
s2 = zeros(37,1);
% r1 = zeros(2,1); r2 = zeros(2,1);

tmpobj = zeros(2,4);
SOW = sow_idx(alpha(1));
latidx = alternormal(alpha(2),100,36);


for k = 1:100
    tmp_idx = SOW(k); 
    tinf = inf_pool(tmp_idx,:);
%     tlat = lat_pool(tmp_idx,:);
    tlat = latidx(k,:);
%     pause = 1;
    s1(:) = 0; s1(1) = 120e8; 
    s2(:) = 0; s2(1) = 25e8;
    
    tote = 0;
    mine = 1e8;
    wd = 0;
    midx = 1;
    tmpmr = 0;
    for i = 1:36
        tmpinf = exp(tinf(i)*std_mm(i) + inf_mm(i));
%         tmplat = b1 + b2*tmpinf + alpha(2)*tlat(i)/4*(db1 + db2*tmpinf);
        tmplat = b1 + b2*tmpinf + tlat(i)/4*(db1 + db2*tmpinf);
        tmpinf = tmpinf*10*24*3600;
        tmplat = tmplat*10*24*3600;
        
        if type == 0
            xvar = [(s1(i)-smin1)/(smax1-smin1) (s2(i)-smin2)/(smax2-smin2) mod(i-1,36)/36 (tinf(i)+2)/4 (tlat(i)+2)/4]; % ref0 xvar
        else
            xvar = [(s1(i)-smin1)/(smax1-smin1) (s2(i)-smin2)/(smax2-smin2) mod(i-1,36)/36 (tinf(i)>0)*0.5+0.25 (tlat(i)>0)*0.5+0.25]; % ref0 xvar
        end
        
        tmpd1 = RBF(xvar,Cp,Rp,Wp1)*(rmax1(i)-rmin1(i)) + rmin1(i);
%         dd1(i) = RBF(xvar,Cp,Rp,Wp1);
        tmpd1 = tmpd1*10*24*3600;
        tmpa1 = max([0 s1(i)+tmpinf-smax1-tmpd1]); 
        tmpr1 = min([tmpa1+tmpd1 s1(i)+tmpinf-smin1]);
        s1(i+1) = s1(i) + tmpinf - tmpr1;
%         r1(tmark) = tmpr1/10/24/3600;
        
        tmpd2 = RBF(xvar,Cp,Rp,Wp2)*(rmax2(i)-rmin2(i)) + rmin2(i);
%         dd2(i) = RBF(xvar,Cp,Rp,Wp2);
        tmpd2 = tmpd2*10*24*3600;
        tmpa2 = max([0 s2(i)+tmplat+tmpr1-smax2-tmpd2]);
        tmpr2 = min([tmpd2+tmpa2 s2(i)+tmplat+tmpr1-smin2]);
        s2(i+1) = s2(i) + tmpr1 + tmplat - tmpr2;
%         r2(tmark) = tmpr2/10/24/3600;
        
        tmpe = 0;
        tmpe = tmpe + 8.5*1e3*min([tmpr1 tmpd1])*(inhead1(s1(i))-outhead1(tmpr1));
        tmpe = tmpe + 8.5*1e3*min([tmpr2 tmpd2])*(inhead2(s2(i))-outhead2(tmpr2));
        for kk = 1:10
            tmpe = tmpe + 8.5*1e3*kkhead(kk)*(tmpr1 + tmplat*kkcsa(kk));
        end
        tmpe = tmpe/1e8/1e3/3600;
        
        mine = min(tmpe, mine);
        tote = tote + tmpe; 
        
        tmpmr = tmpmr + tmpr2/10/24/3600/3;
        if mod(i,3) == 0 
            wd = wd + 1/36*max(0, 1.2*demand(midx) - tmpmr)/1.2/demand(midx);
            midx = midx + 1;
            tmpmr = 0; 
        end     
%         pause = 1;
    end
    sp = (s1(37) + s2(37))/(smax1 + smax2);
    tmpobj(k,:) = [tote mine wd sp];
end

count = 0;
for k = 1:100
    count = count + (tmpobj(k,1)>obj_base(1))*(tmpobj(k,2)>obj_base(2))*(tmpobj(k,3)<obj_base(3))*(tmpobj(k,4)>obj_base(4));
    pause = 1;
end
% obj = mean(tmpobj); 
% obj = tmpobj;

end














