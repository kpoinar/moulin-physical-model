function T_row = solveTmoulin(Tprev,Tpmp,Tfar,C,dFprev,nx,dx,dt)
%
cp = 2115;  % J/kgK
Lf = 3.35e5; % J/kg
%
%
R = Tprev;
R(1) =  Tpmp;               % Dirichlet BC at crevasse wall
R(end) = Tfar;              % Dirichlet BC at far-field point
% ADD THE SOURCE TERM FROM LATENT HEAT:
R(2) = R(2) + Lf/cp * abs(dFprev) ./ dx(2);
% What I have done here is incorrect and sloppy -- I added a source term at
% the 2nd node, so close to the boundary, but not on the boundary.
% In reality, it is a Stefan problem: overdetermined BC (BC is Dirichlet at
% Tpmp as well as a flux boundary condition with the temperature gradient
% such that it can conduct the latent heat in).
% So, this is sloppy and should be thought out better in the future.
%
%
%
%
%
%
% Solve for temperature
% Make matrix diagonals AD, AU, AL
AD = R*0;           % form a vector of zeros of length 2*N+1
AD(2:end-1) = 1 + 2*dt./(dx(1:end-2).*dx(2:end-1)) * C;
AD(1) = 1;                  % Dirichlet BC at crevasse wall
AD(end) = 1;                % Dirichlet BC at far-field point
%
AU = R*0;
AU(2:end-1) = -dt./(dx(1:end-2).*dx(2:end-1)) * C;
AU(end) = AU(end-1);        % This coef doesn't make it into the matrix
AU(2) = 0;                  % Dirichlet BC at lower bndry of ground
%
AL = R*0; 
% AL(1) doesn't matter, doesn't make it into the matrix.
AL(1:end-2) =  -dt./(dx(1:end-2).*dx(2:end-1))* C;
AL(end-1) = 0;              % Dirichlet BC at ground surface (air)
%
%
% Set up matrix and solve for T:
A = spdiags([AL AD AU],[-1 0 1],nx,nx);
T_row = A\R;
%