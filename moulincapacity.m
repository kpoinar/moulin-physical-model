function [egg_cap, water_cap] = moulincapacity(M, z, hw)
% calculate the capacity of the moulin.
    m_ell    = M.r_minor .* M.r_major .* pi;
    m_cir    = pi .* M.r_minor .^2;
    tot_area = m_ell /2 + m_cir /2;  
    egg_cap = trapz(z, tot_area);
    
    [~,j] = min(abs(z-hw));
    water_cap = trapz(z(1:j), tot_area(1:j));
end