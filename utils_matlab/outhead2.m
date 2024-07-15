function Y = outhead2(rr)
    rr = rr/10/3600/24;
    Y = -0.0000002574*rr.^2 + 0.0022264693*rr + 1619.5630307963;
end
