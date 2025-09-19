function A = MatrixA(PRM,V,EPos,Pplus)

A = zeros(PRM.N+3,PRM.N+3);

% BC at x=0
A(1,1) = -PRM.hbar;
A(1,2) = 2*1i*PRM.p*PRM.deltax;
A(1,3) = PRM.hbar;

% BC at x=L
A(PRM.N+3,PRM.N+1) = -PRM.hbar;
A(PRM.N+3,PRM.N+2) = -2*1i*Pplus*PRM.deltax;
A(PRM.N+3,PRM.N+3) = PRM.hbar;

% SE
for i = 2:(PRM.N+2)
    A(i,i-1) = -(PRM.hbar^2)/(2*PRM.mstar);
    A(i,i) = PRM.hbar^2/PRM.mstar-(PRM.q*V(i-1)+EPos)*(PRM.deltax)^2;
    A(i,i+1) = -(PRM.hbar^2)/(2*PRM.mstar);
end

end