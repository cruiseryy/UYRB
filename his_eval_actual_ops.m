clear;clc;close all

addpath('utils_matlab');
addpath('data_matlab');

% load parameters
load('paras.mat')

% read historic 10-day flow data 2001-2009 
load('flow_data_new.mat');

% initialize the LYX and LJX storage states s1, s2
s1 = []; s2 = []; 
s1(1) = 120e8; s2(1) = 25e8;

tmark = 1;
tmpmr = 0; 
water_deficit = zeros(9,1);
hydro_power = [];
storage_penalty = [];
for y = 1:9 
    for i = 1:36 
        
        inf = lyx_in(tmark);
        r1 = lyx_out(tmark);
        lat = ljx_in(tmark) - lyx_out(tmark);
        r2 = ljx_out(tmark);
        
        inf = inf*10*24*3600;
        lat = lat*10*24*3600;
        r1 = r1*10*24*3600;
        r2 = r2*10*24*3600;
        
        s1(tmark+1) = s1(tmark) + inf - r1; 
        s2(tmark+1) = s2(tmark) + lat + r1 - r2; 
        
        % estimate total hydropower production of the 10-day period
        tmp_hp = 0;
        tmp_hp = tmp_hp + 8.5*1e3*r1*(inhead1(s1(tmark)) - outhead1(r1));
        tmp_hp = tmp_hp + 8.5*1e3*r2*(inhead2(s2(tmark)) - outhead2(r2));
        for k = 1:10
            tmp_hp = tmp_hp + 8.5*1e3*kkhead(k)*(r1 + lat*kkcsa(k));
        end
        tmp_hp = tmp_hp/1e8/1e3/3600;
        hydro_power(tmark) = tmp_hp;
        
        % estimate the water deficit but at a monthly scale compatible with the demands
        tmpmr = tmpmr + r2/10/24/3600/3;
        if mod(i, 3) == 0
            midx = floor((i-0.1)/3) + 1;
            water_deficit(y) = water_deficit(y) + 1/36 * max(0, 1.2*demand(midx) - tmpmr) / 1.2 / demand(midx);
            tmpmr = 0;
        end
        
        tmark = tmark + 1;
    end
    
    storage_penalty(y) = (s1(tmark) + s2(tmark))/(smax1+smax2);
end

tot_hp = [];
min_hp = [];
for y = 1:9
    tmp_hp = hydro_power((y-1)*36+1:36*y);
    tot_hp(y) = sum(tmp_hp);
    min_hp(y) = min(tmp_hp);
end

his_obj = [tot_hp' min_hp' water_deficit storage_penalty'];