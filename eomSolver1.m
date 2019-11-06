function dydt = eomSolver1(t,y,ub,c,n0, xi0, xi1, xi2, xi3)

    %this function is continuously called at each step
    delta = xi0+xi1*t+xi2*t^2+xi3*t^3;
    %disp(delta);
    ratio = c*n0/(c-n0.*t);

    dydt = zeros(4,1);
    dydt(1) = -1*(ub-y(3)*y(2)^2)/(y(3)^2)+ratio.*sin(delta);
    dydt(2) = -1*y(1)*y(2)/y(3)+ratio.*cos(delta);
    dydt(3) = y(1);
    dydt(4) = y(2)/y(3);
end