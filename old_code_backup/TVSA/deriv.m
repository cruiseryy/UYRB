function d = deriv(xvar, C, R, W, k)

% 5 RBFs and each takes 5 inputs
na = 5;
nb = 5;

% initialize the derivative of Decision to the kth input
d = 0;

% get the centers, radius & input of the kth input
tC = C(k:na:end);
tR = R(k:na:end);
tx = xvar(k);

for i = 1:na
    tmpd = 0;
%     interate over all inputs to compute the exponent part
    for j = 1:nb
        if abs(R(nb*(i-1)+j)) >= 1e-3
            tmpd = tmpd - (xvar(j)-C(nb*(i-1)+j))^2/R(nb*(i-1)+j)^2;
        end
    end
%     multiply by the derivate part for the kth input
    if abs(tR(i)) >= 1e-3
        d = d + W(i)*(-2)*(tx-tC(i))/tR(i)^2*exp(tmpd);
    end
end

end