% Deformation of moulin due to ice deformation via Glen's Flow Law
%
function dG = deformGlen(H, T, alpha, z, n, dt, C)
%
% Glen's Flow Law integrated through the ice column
%
% Pressure:
P = C.rhoi*C.g*(H-z);
A = flowlawexponent(mean(T,2),P,C) ;
%
% u_defm = 2 * (rho * g * alpha).^n * integrate_bed^sfc ( A(z) * (H-z)^n) dz
%
udefm = abs(2 * (C.rhoi * C.g * alpha).^n  *  cumtrapz(z, A.*(H-z).^n));
%
% Translate to deformation
dG = udefm * dt;
