clear variables
close all
clc

% Numerical solution of the single particle Schrödinger equation
% subject to a RTD-type potential solved by a finite difference scheme
%
% Last update: 19 September 2025
%
% Authors: Giulia Elena Aliffi, University of Catania, Italy,
%               e-mail: giuliaelena.aliffi@phd.unict.it
%          Giovanni Nastasi, University of Enna "Kore", Italy,
%               e-mail: giovanni.nastasi@unikore.it
%          Vittorio Romano, University of Catania, Italy,
%               e-mail: vittorio.romano@unict.it
%
% References:
% [1] A. Arnold, Mathematical concepts of open quantum boundary conditions,
%     Transp. Theory Stat. Phys., 30, 561–584 (2001), https://doi.org/10.1081/TT-100105939
% [2] A. Jüngel, Transport Equations for Semiconductors, Springer, Berlin (2009)
% [3] G.E. Aliffi, G. Nastasi, V. Romano, Ballistic electron transport 
%     described by a fourth-order Schrödinger equation, 
%     Z. Angew. Math. Phys., 76:155 (2025), https://doi.org/10.1007/s00033-025-02535-5
%
% MIT License
% Copyright (c) 2025 Giulia Elena Aliffi, Giovanni Nastasi and Vittorio Romano


PRM = set_parameters;
V = potential(PRM);

figure
plot(PRM.x,V)
axis([0, PRM.L, min(V), max(V)])
xlabel('length (nm)')
ylabel('electrostatic potential (V)')
title('RTD-type potential')

[RealPsi,ImagPsi,Psi] = PsiSinglep(PRM,V);
module2df = sqrt(RealPsi.^2 + ImagPsi.^2);

figure
plot(PRM.x,module2df(2:PRM.N+2))
axis([0, PRM.L, min(module2df(2:PRM.N+2)), max(module2df(2:PRM.N+2))])
xlabel('length [nm]');
ylabel('\Psi');
title('Module of the wavefunction')

J2df = CurrentSinglep(PRM,Psi);

figure
plot(PRM.x,J2df)
axis([0, PRM.L, min(J2df), max(J2df)])
xlabel('length [nm]');
ylabel('probability current density [fs^-^1 nm^-^2]');
title('Single electron current')