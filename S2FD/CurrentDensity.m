function J=CurrentDensity(PRM,V)

nn =[0:(PRM.n-2)]'; % Crea il vettore in base all'input
alpha= zeros(1,PRM.n); % rispetto al libro sta scambiando alpha e beta
beta = (nn+1)./sqrt(4*(nn+1).^2 -1) ; % Determina le beta 
T = diag(alpha) + diag(beta,1) + diag(beta,-1); % costruisce T
[autovettori,autovalori]=eig(T); % Determina autovalori e autovettori della matrice T
nodi=diag(autovalori) ; % nodi in [-1,1]
z_i=autovettori(1,:)'; % Seleziona in un vettore colonna le prime componenti degli autovettori normalizzati
w=(z_i.^2)*2; % Determina i pesi in ab
c=(PRM.b-PRM.r)/2;
d=(PRM.b+PRM.r)/2;
pp=nodi.*c + d; % determina i nodi in ab

EP=zeros(length(pp),1);

for v=1:length(pp)

    EP(v)=Energy(PRM,pp(v));

end

R=zeros(PRM.n,PRM.N+3); % con punti ghost
I=zeros(PRM.n,PRM.N+3);

for i=1:length(pp)

    [RealPsi,ImagPsi,~]=PsiSinglep(PRM,V,pp(i));
    R(i,:)=RealPsi';
    I(i,:)=ImagPsi';

end

CONJ=R-1i*I;
PSI=R+1i*I;

% Definisco la matrice delle derivate nei punti griglia

PSIPRIME=zeros(PRM.n,PRM.N+1); % ! non ho punti ghost

for z=1:PRM.n
    for j=1:(PRM.N+1)
        
        PSIPRIME(z,j)=(PSI(z,j+2)-PSI(z,j))/(2*PRM.deltax);
        
    end
end

% Calcolo la parte immaginaria di psicong*psiprime nei punti griglia
% interni

IMPSI=imag(CONJ(1:PRM.n,2:(PRM.N+2)).*PSIPRIME);

gw=log(1+exp((PRM.EF-EP)/PRM.KBT)).*w; 
f=gw*ones(1,PRM.N+1).*IMPSI;
quadratura=c*sum(f,1);
% Current=PRM.q*PRM.KBT*quadratura/(2*pi^2*PRM.hbar^2);
% plot(PRM.x,Current);  % Se la differenza di potenziale è nulla, mi aspetto che venga 0

electron=1.602e-19; % verifica
J=-electron*PRM.KBT*quadratura/(2*pi^2*PRM.hbar^2)*1.e12; % unità ampere/micron^2 [A/um^2]

    
end