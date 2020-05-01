%% Initial values and coordinate grid
clear all ;
x_s = -1.62 ;
y_s = 0 ; 
z_s = 3.06 ; 

x_m = 0.38 ; 
y_m = 0 ; 
z_m = 0.46 ; 

d = 0.4 ;

dx = 0.05 ; dy = 0.05 ; dz = 0.05 ;
x = [-0.3:dx:0.3] ; y = [-0.3:dy:0.3] ; z = [0.35:dz:3.15] ;
sz_x = length(x) ; sz_y = length(y) ; sz_z = length(z) ;

[X,Y,Z] = meshgrid(x,y,z) ;

c0 = 343 ;
f_array = 100:100:3000 ;

V = 20 ;

step = 0.02 ;

%% Integrating  
X = X(:) ; Y = Y(:) ; Z = Z(:) ;
cub = length(X) ;

dist_x = X * ones(1,cub) - ones(cub,1)*X' ;
dist_y = Y * ones(1,cub) - ones(cub,1)*Y' ;
dist_z = Z * ones(1,cub) - ones(cub,1)*Z' ;


dist = sqrt(dist_x.^2 + dist_y.^2 + dist_z.^2) ;
clear dist_x dist_y dist_z

r2_min_r1 = sqrt((x_s-x_m)^2 + (y_s-y_m)^2 + (z_s-z_m)^2) ;
r2_min_r  = sqrt((X-x_m).^2 + (Y-y_m).^2 + (Z-z_m).^2) ;
r1_min_r  = sqrt((X-x_s).^2 + (Y-y_s).^2 + (Z-z_s).^2) ;    
r = sqrt(X.^2 + Y.^2 + Z.^2) ;

%% Calculating Mach number
Mxy = Mach(X, Y, Z, V, d, dz) ;
M = reshape(Mxy, sz_x, sz_y, sz_z) ;

Mxy_mat = diag(Mxy) ;
p_sum = [] ;
clear X Y Z

%% Buiding 3D gradient and laplacian matrices

D_3d = -diag(ones(1, cub), 0)/dx - diag(ones(1, cub), 0)/dy - ...
    diag(ones(1, cub), 0)/dz + diag(ones(1, cub-1), 1)/dx + ...
    diag(ones(1, cub-sz_x), sz_x)/dy + diag(ones(1, cub-sz_x*sz_y), sz_x*sz_y)/dz ;
D_3d((1:sz_x*sz_y)*sz_y, (1:sz_x*sz_y)*sz_y + 1) = 0 ;
D_3d = D_3d(1:cub, 1:cub) ;


D2_3d = -2*diag(ones(1, cub), 0)/dx^2 - 2*diag(ones(1, cub), 0)/dy^2 -...
    2*diag(ones(1, cub), 0)/dz^2 ;
D2_3d = D2_3d + diag(ones(1, cub-1), 1)/dx^2 + ...
    diag(ones(1, cub-1), -1)/dx^2 + diag(ones(1, cub-sz_x), -sz_x)/dy^2;
D2_3d = D2_3d + diag(ones(1, cub-sz_x), sz_x)/dy^2 + diag(ones(1, cub-sz_x*sz_y), sz_x*sz_y)/dz^2 +...
    diag(ones(1, cub-sz_x*sz_y), -sz_x*sz_y)/dz^2 ;

D2_3d((1:sz_x*sz_y)*sz_y + 1, (1:sz_x*sz_y)*sz_y) = 0 ;
D2_3d((1:sz_x*sz_y)*sz_y, (1:sz_x*sz_y)*sz_y + 1) = 0 ;
D2_3d = D_3d(1:cub, 1:cub) ;

    
%% Solving the integral equation
for f = f_array
    f
    omega = 2*pi*f ;
    k = omega/c0 ;
    r_array = (step/10 : step/10 : 4) ;
    G = 1/(4*pi)* exp(-1i*k*r_array)./r_array ;
    r_array = [0, r_array] ;
    G = [1, G] ;
    G_r = interp1(r_array, G, dist) ;
    
    psi0 = 1/(8*pi^3*c0^2)*interp1(r_array, G, r2_min_r1)* ones(cub,1) ;
    coeff = (-2i*Mxy_mat*k.*D_3d + Mxy_mat.^2.*D2_3d) ;
    mat_coef = - coeff .* G_r + eye(cub); 
    clear coeff G G_r
    p = (mat_coef) \ psi0 /(dx*dy*dz);
    p_sum = [p_sum, sum(p)] ;        
end %end for each frequency

figure; 
hold all;
plot(real(p_sum)) ;
plot(imag(p_sum)) ;
plot(abs(p_sum)) ;
save('data_calc_0') ;
    