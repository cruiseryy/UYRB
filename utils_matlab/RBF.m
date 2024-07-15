function d = RBF(xvar, C, R, W)
d = 0;
na = 5;
nb = 5;
for i = 1:na
    tmpd = 0;
    for j = 1:nb
        if R(nb*(i-1)+j) ~= 0
            tmpd = tmpd - (xvar(j)-C(nb*(i-1)+j))^2 / R(nb*(i-1)+j)^2;
        end
    end
    d = d + W(i)*exp(tmpd);
end

end 