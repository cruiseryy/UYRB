% clear;clc;close all
% clear;clc
% rng(0)

function idx = sow_idx(b)

% threshold
alpha = 0.25;

% sample size 
n = 100; 

% scaling factor 
% b = 5; 

inf_mm = [5.07559,4.98739,5.02828,5.03286,5.09163,5.17904,5.22714,5.36965,5.43714,5.57579,5.64422,5.87030,5.90242,5.94810,6.14867,6.38404,6.54943,6.78031,6.87444,6.90823,6.80754,6.72528,6.76357,6.78094,6.86762,6.78863,6.83915,6.88035,6.72561,6.49341,6.23301,6.01686,5.73047,5.51399,5.35356,5.21524];
std_mm = [0.34692,0.30352,0.23062,0.20955,0.22309,0.20396,0.16340,0.16398,0.19826,0.16088,0.21378,0.30881,0.37468,0.32197,0.26952,0.27893,0.32660,0.39741,0.31701,0.45362,0.45336,0.42666,0.41728,0.47373,0.50394,0.41208,0.35528,0.46297,0.45632,0.38720,0.32151,0.29521,0.28154,0.24770,0.26209,0.22290];

inf_pool = importdata('../streamflow_generator/inflow.txt');
% tinf = inf_pool(k,:);
% tmpinf = exp(tinf(i)*std_mm(i) + inf_mm(i));

total_inf = exp(inf_pool.*repmat(std_mm,size(inf_pool,1),1) + repmat(inf_mm,size(inf_pool,1),1));
total_ainf = mean(total_inf,2);

% lower and upper total annual runoff for the threshold alpha
lowth = quantile(total_ainf,alpha);
lowidx = find(total_ainf<=lowth);
uppth = quantile(total_ainf,1-alpha);
uppidx = find(total_ainf>=uppth);
mididx = find((total_ainf>lowth).*(total_ainf<uppth)==1);

sample_ainf = zeros(n,1);
idx = zeros(n,1);
for k = 1:n
    u = rand();
    tmp_case = 1*(u<=b/(2+2*b)) + 2*(u>b/(2+2*b)&&u<=2*b/(2+2*b)) + 3*(u>2*b/(2+2*b));
    switch tmp_case
        case 1
            tmp_infidx = lowidx(randi(5000));
        case 2
            tmp_infidx = uppidx(randi(5000));
        case 3
            tmp_infidx = mididx(randi(10000));
    end
    sample_ainf(k) = total_ainf(tmp_infidx);
    idx(k) = tmp_infidx;
end

end

% figure()
% plot(mean(total_inf,1))

















