function idx = sow_idx(af, n)

% a threshold for defining 'extreme' dry conditions
% i use the 10th percentile of total flow here
alpha = 0.1;

load('paras.mat', 'inf_mm')
load('paras.mat', 'std_mm')

inf_pool = importdata('inflow.txt');

total_inf = exp(inf_pool.*repmat(std_mm,size(inf_pool,1),1) + repmat(inf_mm,size(inf_pool,1),1));
total_ainf = sum(total_inf,2);

dry_th = quantile(total_ainf, alpha);
dry_idx = find(total_ainf <= dry_th);
normal_idx = find(total_ainf > dry_th);

idx = zeros(n,1);
for k = 1:n
    u = rand();
    if u <= af*alpha
        idx(k) = datasample(dry_idx, 1);
    else
        idx(k) = datasample(normal_idx, 1);
    end
end

end

















