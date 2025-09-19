function Jp = CurrentSinglep(PRM,Psi)

ConjPsi=conj(Psi);
Jp=zeros(1,PRM.N+1);

for j=0:PRM.N
    A=PRM.hbar/(PRM.mstar*2*PRM.deltax);
    Jp(j+1)=A*imag(ConjPsi(j+2)*(Psi(j+3)-Psi(j+1)));
end

end