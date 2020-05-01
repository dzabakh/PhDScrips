% Building Mach number M(x,y,z)
function [Mxy] = Mach(Xgrid, Ygrid, Zgrid, U, r_nozzle, dz)
    a = 0.07 ; % turbulence parameter
    c0 = 343 ;
    
    tan_beta = 3.4 * a ; % angle of divergence
    x_per = 5.96*r_nozzle/a ; % flux pole to transitional cross-section
    L0 = 0.29*r_nozzle/a ;    % flux pole to initial cross-section

    Ygrid = Ygrid + L0 ;

    z2 = Ygrid(find(Ygrid >= L0 & Ygrid <= x_per)) ;
    z3 = Ygrid(find(Ygrid > x_per)) ;
    z = [z2; z3] ;
    
    r_yadra = z*0 ;
    r_yadra(find(Ygrid >= L0 & Ygrid <= x_per)) = ...
    r_nozzle - r_nozzle/(x_per - L0).*(z2-L0) ;
  
    r_stream = Ygrid(find(Ygrid >= L0))*tan_beta + 0.01 ;
    
%% x > x_per
%%% amplitude
    L = Ygrid - L0;
    um_div_u0 = 0.96*r_nozzle/a./(L + 0.29*r_nozzle/a) ;
    um_div_u0(z < L0) = 1 ;
%% L0 < x < x_per
    um_div_u0(L0+dz < Ygrid < 2.908) = 1/2 .* ...
        (1 - tanh((Ygrid - 3.47)/(r_nozzle * 2))) ;
    
%     um_div_u0 = 1 ;
%     Mxy = um_div_u0*U/c0 * 0.5 .* ...
%     (1 - tanh((sqrt(Xgrid.^2 + Zgrid.^2)-r_nozzle)/r_nozzle*10));    
    
    Mxy = um_div_u0*U/c0 * 0.5 .* ...
    (1 - tanh((sqrt(Xgrid.^2 + Zgrid.^2) - r_yadra)./...
    (r_stream-r_yadra)/0.2));

%     load('stream_data.mat')
%     
%     Ms = reshape(Mxy, 70,70,70);
%     figure; hold all; 
%     plot((X+305)/1000+0.4, Ux_0/343);
%     plot((1:70)*0.02,Ms(:,1,35)); 
%     figure; hold all; 
%     plot((X+305)/1000+0.4, Ux_30/343);
%     plot((1:70)*0.02,Ms(:,15,35)); 
%     figure; hold all;
%     plot((1:70)*0.02,Ms(:,30,35)); 
%     plot((X+305)/1000+0.4, Ux_60/343);
 
end