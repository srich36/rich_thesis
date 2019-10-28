function dydt = odefun(t,y,ub,c,n0, xi0, xi1, xi2, xi3)

    deltas = deltaT1(t, xi0, xi1, xi2, xi3);

    dydt = zeros(4,1);
    dydt(1) = (ub-y(3)*y(2)^2)/(y(3)^2)+thArc1ThrustRatio(c, n0, t).*sin(deltas);
    dydt(2) = y(1)*y(2)/y(3)+thArc1ThrustRatio(c, n0, t).*cos(deltas);
    dydt(3) = y(1);
    dydt(4) = y(2)/y(3);
end