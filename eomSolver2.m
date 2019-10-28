function dydt = eomSolver2(t,y,ub,c,n0, t1, t2, v0, v1, v2, v3)

    deltas = deltaT2(t, t2, v0, v1, v2, v3);

    dydt = zeros(4,1);
    dydt(1) = (ub-y(3)*y(2)^2)/(y(3)^2)+thArc1ThrustRatio(c, n0, t1, t, t2).*sin(deltas);
    dydt(2) = y(1)*y(2)/y(3)+thArc1ThrustRatio(c, n0, t1, t, t2).*cos(deltas);
    dydt(3) = y(1);
    dydt(4) = y(2)/y(3);

end