function [RealPsi,ImagPsi,Psi]=PsiSinglep(PRM,V)

p = PRM.p;
if p>0
    if (p^2+2*PRM.mstar*PRM.q*(PRM.VL-PRM.V0))>0
        Pplus = sqrt(p^2+2*PRM.mstar*PRM.q*(PRM.VL-PRM.V0));
        EPos = p^2/(2*PRM.mstar)-PRM.q*PRM.V0;
        A = MatrixA(PRM,V,EPos,Pplus);
        b = zeros(PRM.N+3,1);
        b(1) = 4*1i*p*PRM.deltax;
        Psi = A\b;
        RealPsi = real(Psi);
        ImagPsi = imag(Psi);
    else
        Pplus = 1i*sqrt(-p^2-2*PRM.mstar*PRM.q*(PRM.VL-PRM.V0));
        EPos = p^2/(2*PRM.mstar)-PRM.q*PRM.V0;
        A = MatrixA(PRM,V,EPos,Pplus);
        b = zeros(PRM.N+3,1);
        b(1) = 4*1i*p*PRM.deltax;
        Psi = A\b;
        RealPsi = real(Psi);
        ImagPsi = imag(Psi);
    end
else
    if (p^2-2*PRM.mstar*PRM.q*(PRM.VL-PRM.V0))>0
        Pminus = sqrt(p^2-2*PRM.mstar*PRM.q*(PRM.VL-PRM.V0));
        ENeg = p^2/(2*PRM.mstar)-PRM.q*PRM.VL;
        C = MatrixC(PRM,V,ENeg,Pminus);
        d = zeros(PRM.N+3,1);
        d(PRM.N+3) = 4*1i*p*PRM.deltax;
        Psi = C\d;
        RealPsi = real(Psi);
        ImagPsi = imag(Psi);
    else
        Pminus = 1i*sqrt(-p^2+2*PRM.mstar*PRM.q*(PRM.VL-PRM.V0));
        ENeg = p^2/(2*PRM.mstar)-PRM.q*PRM.VL;
        C = MatrixC(PRM,V,ENeg,Pminus);
        d = zeros(PRM.N+3,1);
        d(PRM.N+3) = 4*1i*p*PRM.deltax;
        Psi = C\d;
        RealPsi = real(Psi);
        ImagPsi = imag(Psi);
    end
end

end