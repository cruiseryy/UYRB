function Y = inhead2(ss)
    ss = ss/100000000;
    Y = -0.0217686085*ss.^2 + 2.1651210663*ss + 1681.8568952269;
end