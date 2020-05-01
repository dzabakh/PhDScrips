% Building Mach number M(x,y,z)

function [Mxy] = Mach(Xgrid, Ygrid, Zgrid, U, d, dz)
%     d = 0.4 ;
%     U = 50 ;
    a = 0.07 ; % turbulence parameter
    c0 = 343 ;
%     dz = 0.05 ;
%     x = [-0.5:0.05:0.5] ;
%     y = [-0.5:0.05:0.5] ;
%     z = [dz:dz:3.5] ;
    
%     sz_x = length(x) ;
%     sz_y = length(y) ;
%     sz_z = length(z) ;
    
%     [Xgrid, Ygrid, Zgrid] = meshgrid(x,y,z) ;
% 
%     Xgrid = Xgrid(:) ;
%     Ygrid = Ygrid(:) ;
%     Zgrid = Zgrid(:) ;
    
    tan_beta = 3.4 * a ; % angle of divergence
    x_per = 0.96*d/2/a ; % flux pole to transitional cross-section
    L0 = 0.29*d/2/a ;    % flux pole to initial cross-section


    Zgrid = Zgrid + L0 ;


    z1 = Zgrid(find(Zgrid <= L0)) ;
    z2 = Zgrid(find(Zgrid > L0 & Zgrid <= x_per)) ;
    z3 = Zgrid(find(Zgrid > x_per)) ;


    r1 = z1 * tan_beta ;
    r2 = d/2 - d/2/(x_per - L0).*(z2-L0) ;
    r3 = z3 * 0 ;
    r_yadra = [r2; r3] ;
    z = [z2; z3] ;
    figure;
    plot(z,r_yadra);
    hold all

    r_stream = z*tan_beta ;

    plot(Zgrid, r_stream);


%% x > x_per

%%% amplitude
    r0 = d/2 ;
    L = Zgrid - L0;
    um_div_u0 = 0.96*r0/a./(L + 0.29*r0/a) ;
    um_div_u0(z < L0) = 1 ;

%% L0 < x < x_per
    um_div_u0(L0+dz < Zgrid < 2.908) = 1/2 .* ...
        (1 - tanh((Zgrid - 3.47)/d)) ;


%     end    

% [Ymesh, Xmesh] = meshgrid(Y, z) ;
% figure;



    figure; 
    plot(z, um_div_u0);
    Mxy = um_div_u0*U/c0 * 0.5 .* ...
    (1 - tanh((sqrt(Xgrid.^2 + Ygrid.^2) - r_yadra)./...
    (r_stream-r_yadra)));


end
