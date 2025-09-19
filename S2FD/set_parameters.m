function PRM = set_parameters 

% Physical parameters
PRM.hbar = 0.6582;              % Reduced Planck constant [eV fs]
PRM.q = 1;                      % Elementary charge [eV/V]
PRM.mstar = 0.380937;           % Effective mass of GaAs [eV fs^2 / nm^2]

% Particle parameters
PRM.p = 0.2;                    % Momentum of the particle  [eV fs / nm]

% Device parameters
PRM.L = 135;                    % Length [nm] 
PRM.V0 = 0;                     % Applied voltage at x=0
PRM.VL = 0.1;                   % Applied voltage at x=L 
PRM.a1 = 50;                    % Length of the 1st junction [nm] 
PRM.a2 = 60;                    % Length of the 2nd junction [nm] 
PRM.a3 = 65;                    % Length of the 3rd junction [nm] 
PRM.a4 = 70;                    % Length of the 4th junction [nm] 
PRM.a5 = 75;                    % Length of the 5th junction [nm] 
PRM.a6 = 85;                    % Length of the 6th junction [nm] 
PRM.vbar = -0.3;                % Barrier potential [V]

% Numerical parameters
PRM.N = 5000;                   % Number of cells
PRM.deltax = PRM.L/PRM.N;       % Mesh size
PRM.x = 0:PRM.deltax:PRM.L;     % Mesh points

end