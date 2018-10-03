function [x,dx,nx] = setupx(dt,chebx,xmax)

%Constants
ki =    2.1; % J/mKs
rhoi = 910; % kg/m3
cp = 2115;  % J/kgK
C = ki /(rhoi * cp); 
% The diffusion lengthscale:
dx = sqrt(C * dt);   % meters
switch chebx
    case 0
        % uniform resolution in x
        x = 0:dx:xmax;
    case 1
        % concentrate x resolution at sidewalls
        x = linspace(0,pi/2,round(xmax/dx/10));
        x = (xmax * (1-sin(x)));
    otherwise 
        error('chebx: enter 0 or 1')
end
%
nx = length(x);
dx = diff(x);
%dx = [dx(1) dx];
dx = [dx dx(end)];