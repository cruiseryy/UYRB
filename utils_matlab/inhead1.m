function Y = inhead1(ss)
    ss = ss/100000000;
    Y = -0.0007252553*ss.^2 + 0.5690388196*ss + 2502.7122406983;
end