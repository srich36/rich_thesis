function delta = deltaT2(t, t2, v0, v1, v2, v3)
    delta = v0+v1*(t-t2)+v2*(t-t2)^2+v3*(t-t2)^3;
end