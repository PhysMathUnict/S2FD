function C = MatrixC(PRM,V,ENeg,Pminus)

C = zeros(PRM.N+3,PRM.N+3);

% BC at x=0
C(1,1) = -PRM.hbar;
C(1,2) = 2*1i*Pminus*PRM.deltax;
C(1,3) = PRM.hbar;

% BC at x=L
C(PRM.N+3,PRM.N+1) = -PRM.hbar;
C(PRM.N+3,PRM.N+2) = 2*1i*PRM.p*PRM.deltax;
C(PRM.N+3,PRM.N+3) = PRM.hbar;

% SE
for i = 2:(PRM.N+2)
    C(i,i-1) = -(PRM.hbar^2)/(2*PRM.mstar);
    C(i,i) = PRM.hbar^2/PRM.mstar-(PRM.q*V(i-1)+ENeg)*(PRM.deltax)^2;
    C(i,i+1) = -(PRM.hbar^2)/(2*PRM.mstar);
end

end