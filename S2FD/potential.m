function V = potential(PRM)

V = zeros(1,PRM.N+1);
m1=(PRM.V0-PRM.VL)/(PRM.a1-PRM.a6);
q1=(PRM.a1*PRM.VL-PRM.a6*PRM.V0)/(PRM.a1-PRM.a6);

for z = 1:(PRM.N+1)
    if  PRM.x(z)>=0 && PRM.x(z)<PRM.a1
        V(z) = PRM.V0;
    end
    if PRM.x(z)>=PRM.a1 && PRM.x(z)<PRM.a2
        V(z) = m1*PRM.x(z) + q1;
    end
    if PRM.x(z)>=PRM.a2 && PRM.x(z)<PRM.a3
        V(z) = m1*PRM.x(z) + q1 + PRM.vbar;
    end
    if PRM.x(z)>=PRM.a3 && PRM.x(z)<PRM.a4
        V(z) = m1*PRM.x(z) + q1;
    end
    if PRM.x(z)>=PRM.a4 && PRM.x(z)<PRM.a5
        V(z) = m1*PRM.x(z) + q1 + PRM.vbar;
    end
    if PRM.x(z)>=PRM.a5 && PRM.x(z)<PRM.a6
        V(z) = m1*PRM.x(z) + q1;
    end
    if PRM.x(z)>=PRM.a6 && PRM.x(z)<=PRM.L
        V(z) = PRM.VL;
    end
end

end