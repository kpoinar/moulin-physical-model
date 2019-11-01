function Tz = importTz(Tdatatype,z)
H = max(z);
switch Tdatatype
    case 'Luthi'
        % Use Site D borehole
        load fielddata/Luthidata.mat
        % Put it on the 0 to 1 / 0 to H z grid
        Luthi.zetanew = flipud(linspace(Luthi.zeta(1),Luthi.zeta(end),100));
        % ^ zeta here does not go quite all the way 0 to 1
        Luthi.Tsmooth = fastsmooth(interp1(Luthi.zeta,Luthi.T,Luthi.zetanew),10,3,1);
        Tz = interp1(linspace(1,0,100),Luthi.Tsmooth,z/H);
        %
    case 'Temperate'
        % Temperate ice (no refreezing)
        Tz = G.To*ones(size(z));
    case 'Cool'
        Tmin = -5;
        Tz = linspace(0,Tmin,numel(z))' + G.To;
    case 'Cold'
        Tmin = -3;
        Tz = fastsmooth([linspace(0,Tmin,numel(z)/2) linspace(Tmin,Tmin/2,numel(z)/2+1)]',40,3,1) + G.To;
    case 'HarrS2A'
        load fielddata/Harrington_temps_2015.mat
        harr15.S2_A_zeta = (harr15.S2_A_depth_m - min(harr15.S2_A_depth_m)) / (-min(harr15.S2_A_depth_m));
        harr15.S2_A_zeta(1) = 1;
        harr15.Tsmooth = fastsmooth(interp1(harr15.S2_A_zeta,harr15.S2_A_temp_C+273.15,linspace(0,1,100)),10,3,1);
        Tz = interp1(linspace(0,1,100),harr15.Tsmooth,z/H);    
    case 'HarrS4C'
        load fielddata/Harrington_temps_2015.mat
        harr15.S4_C_zeta = (harr15.S4_C_depth_m - min(harr15.S4_C_depth_m)) / (-min(harr15.S4_C_depth_m));
        harr15.S4_C_zeta(1) = 1;
        harr15.Tsmooth = fastsmooth(interp1(harr15.S4_C_zeta,harr15.S4_C_temp_C+273.15,linspace(0,1,100)),10,3,1);
        Tz = interp1(linspace(0,1,100),harr15.Tsmooth,z/H);
end