function [L, alpha, Qinrange, Qinbase] = makeicesheetgeom(H)

% a function to calculate the idealized ice sheet profile for a given
% ice thickness
%%
tau = 100e3;
C = makeConstants;

L = (H.^2).*(C.g*C.rhoi)./(2.*tau);

alpha = tau./(C.rhoi.*C.g.*H);

Qinrange = zeros(size(L));
Qinbase = zeros(size(L));


% Interpolate from Lauren's values on the google doc table
Ltable = (20:10:110)*1e3;
Qmeantable = [3.7 3.2 4.5 4.5 4.2 4.6 3.2 3.4 3.4 3.7];
Qrangetable = [5 2.3 3.2 2.1 1.7 1.1 0.8 0.7 0.7 0.8];

% Qinval is actually the diurnal range
Qinrange = interp1(Ltable,Qrangetable,L,'spline');

% Qmultiplier is the DC base
Qinbase = interp1(Ltable,Qmeantable,L,'spline');

%{
for ii=1:length(L)
    H1 = H(ii);
    if H1 <= 400
        Qinval(ii)         = 2;
        Qmultiplier2(ii) = 0;
        
    elseif H1 >400 && H1 <= 500
        Qinval(ii)         = 2;
        Qmultiplier2(ii) = 1;
        
    elseif H1 > 500 && H1 <= 600
        Qinval(ii)         = 3;
        Qmultiplier2(ii) = 1;
        
    elseif H1 > 600 && H1 <= 700
        Qinval(ii)         = 3;
        Qmultiplier2(ii) = 2;
        
    elseif H1 > 700 && H1 <= 800
        Qinval(ii)         = 3;
        Qmultiplier2(ii) = 3.5;
        
    elseif H1 > 800 && H1 <= 900
        Qinval(ii)         = 5;
        Qmultiplier2(ii) = 4;
        
    elseif H1 > 900 && H1 <= 1000
        Qinval(ii)         = 5;
        Qmultiplier2(ii) = 5;
        
    elseif H1 > 1000 && H1 <= 1100
        Qinval(ii)         = 6;
        Qmultiplier2(ii) = 5;
        
        
    elseif H1 > 1100 && H1 <= 1200
        Qinval(ii)         = 7;
        Qmultiplier2(ii) = 3;
        
        
    elseif H1 > 1200 && H1 <= 1300
        Qinval(ii)         = 8;
        Qmultiplier2(ii) = 3;
        
    elseif H1 > 1300 && H1 <= 1400
        Qinval(ii)         = 9;
        Qmultiplier2(ii) = 3;
        
    else
        disp('This combination might not work...')
        Qinval(ii)         = 9;
        Qmultiplier2(ii) = 3;
        
    end
    
end
%}

