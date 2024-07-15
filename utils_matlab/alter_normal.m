function  z = alter_normal(af, nx, ny)

alpha = 0.1;

% sample size 
n = nx*ny; 
z = zeros(n,1);

for k = 1:n
    u = rand();
    if u <= alpha*af
        z(k) = norminv(0.1*rand());
    else
        z(k) = norminv(0.1 + 0.9*rand());
    end
end
z = reshape(z, nx, ny);
end