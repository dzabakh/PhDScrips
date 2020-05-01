% Building Mach number M(x,y,z)

a = 0.07 ; % turbulence parameter
d = 0.6 ;  % nozzle diameter
U = 20 ;
c0 = 343 ;

tan_beta = 3.4 * a ; % angle of divergence
x_per = 0.96*d/2/a ; % flux pole to transitional cross-section
L0 = 0.29*d/2/a ; % flux pole to initial cross-section
dx = 0.05 ;

len_pole = floor(L0/dx) ;
len_trans = floor((x_per - L0)/dx) ;

x1 = 0:dx:L0 ;
x2 = (L0 + dx):dx:x_per  ;
x3 = (x_per + dx):dx:(x_per + dx*50) ;

r1 = x1 * tan_beta ;
r2 = d/2 - d/2/(x_per - L0).*(x2-L0) ;
r3 = x3* 0 ;
r_yadra = [r1, r2, r3] ;
x = [x1, x2, x3] ;
figure;
plot(x,r_yadra);
hold all

r_stream = x*tan_beta ;

plot(x, r_stream);


%% x > x_per

%%% amplitude
r0 = d/2 ;
L = x - L0;
um_div_u0 = 0.96*r0/a./(L + 0.29*r0/a) ;
um_div_u0(x<L0) = 1 ;

%% L0 < x < x_per
um_div_u0(L0+dx<x<2.908) = 1/2 .* (1 - tanh((x - 3.47)/d)) ;

figure; 
plot(x, um_div_u0);
Y = 0:0.05:4 ;
M = [] ;
for i = 1:length(x)
    r_i = r_stream(i) ;
    r_yadra_i = r_yadra(i) ;
    M = [M; um_div_u0(i)*U/c0 * 0.5 .* (1 - tanh((Y - r_yadra_i)/(r_i-r_yadra_i)))] ;
end    

[Ymesh, Xmesh] = meshgrid(Y, x) ;
figure;
surf(Xmesh, Ymesh, M) ;

M2 = fliplr(M) ;
M_all = [M2, M] ;
Ymesh2 = -fliplr(Ymesh) ;
Ymesh_all = [Ymesh2, Ymesh] ;
Xmesh_all = [Xmesh, Xmesh] ;
surf(Xmesh_all, Ymesh_all, M_all) ;
