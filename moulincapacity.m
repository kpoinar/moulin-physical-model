function [egg_cap, water_cap] = moulincapacity(M, z, wet)
% calculate the capacity of the moulin.
    m_ell    = M.r_minor .* M.r_major .* pi;
    m_cir    = pi .* M.r_minor .^2;
    tot_area = m_ell /2 + m_cir /2;  
    egg_cap = trapz(z, tot_area);
    
    if length(find(wet))>1
        water_cap = trapz(z(wet), tot_area(wet));
    else
        water_cap = 0;
    end
end