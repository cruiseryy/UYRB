function y = kge(fit, true)
    r = corr(true, fit);
    alpha = std(fit) / std(true);
    beta = mean(fit) / mean(true);
    y = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2);
end