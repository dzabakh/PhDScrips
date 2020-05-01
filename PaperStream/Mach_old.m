% Building Mach number M(x,y,z)
function [Mxy] = Mach(Xgrid, Ygrid, Zgrid, U, d, dz)
    a = 0.07 ; % turbulence parameter
    c0 = 343 ;
    
    tan_beta = 3.4 * a ; % angle of divergence
    z_per = 0.96*d/2/a ; % flux pole to transitional cross-section
    L0 = 0.29*d/2/a ; % flux pole to initial cross-section

    Zgrid = Zgrid + L0 ;
    
    len_pole = floor(L0/dz) ;
    len_trans = floor((z_per - L0)/dz) ;

%     r_yadra = d/2 ;
    r_yadra = Zgrid * 0;
    r_yadra(find(Zgrid <= L0)) = Zgrid(find(Zgrid <= L0)) * tan_beta ;
    r_yadra(find(Zgrid > L0 & Zgrid <= z_per)) = ...
        d/2 - d/2/(z_per - L0).*(Zgrid(find(Zgrid > L0 & Zgrid <= z_per))-L0) ;

    r_stream = Zgrid*tan_beta ;

    %% x > x_per
    %%% amplitude
    r0 = d/2 ;
    L = Zgrid - L0;
    um_div_u0 = 0.96*r0/a./(L + 0.29*r0/a) ;
    um_div_u0(Zgrid < L0) = 1 ;
    %% L0 < x < x_per
    um_div_u0(L0 + dz < Zgrid < 2.908) = ...
        1/2 .* (1 - tanh((Zgrid - 3.47)/d)) ;

    M = [] ;

    Mxy = um_div_u0*U/c0 * 0.5 .* ...
    (1 - tanh((sqrt(Xgrid.^2 + Ygrid.^2) - r_yadra)./...
    (r_stream-r_yadra)));
end    
