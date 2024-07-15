function  z = alternormal(b, nx, ny)


% sample size 
n = nx*ny; 
z = zeros(n,1);
tmp_case = zeros(n,1);
load('for_alternormal.mat')
for k = 1:n
    u = rand();
    tmp_case(k) = 1*(u<=b/(2+2*b)) + 2*(u>b/(2+2*b)&&u<=2*b/(2+2*b)) + 3*(u>2*b/(2+2*b));
    switch tmp_case(k) 
        case 1
            tmp_infidx = lowidx(randi(numel(lowidx)));
        case 2 
            tmp_infidx = uppidx(randi(numel(uppidx)));
        case 3 
            tmp_infidx = midix(randi(numel(midix)));
    end
    z(k) = enormal(tmp_infidx); 
end
z = reshape(z, nx, ny);
% pause = 1;
end