function ratio = thArc2ThrustRatio(c, n0, t1, t, t2)
    ratio = c*n0/(c-n0.*(t1+t-t2));
end