function Y = outhead1(rr)
    rr = rr/10/3600/24;
    Y = -0.0000002260*rr.^2 + 0.0034508966*rr + 2449.9359789846;
end