%% Initial values and coordinate grid
clear all ;
x_m = 0.2 ;
y_m = 0.2;
z_m = 0.2 ;

x_s = -0.2 ;
y_s = 0.2 ;
z_s = 1 ;

d = 0.1 ;

dx = 0.05 ; dy = 0.05 ; dz = 0.05 ;
x = [-0.5:dx:0.5] ; y = [-0.5:dy:0.5] ; z = [0.1:dz:1.5] ;
sz_x = length(x) ; sz_y = length(y) ; sz_z = length(z) ;

cub = sz_x * sz_y * sz_z ;
[X,Y,Z] = meshgrid(x,y,z) ;

c0 = 343 ;
f_array = 10:10:3500 ;

V = 50 ;

%% Buiding 3D gradient and laplacian matrices

D_3d = -diag(ones(1, cub), 0)/dx - diag(ones(1, cub), 0)/dy - ...
    diag(ones(1, cub), 0)/dz + diag(ones(1, cub-1), 1)/dx + ...
    diag(ones(1, cub-sz_x), sz_x)/dy + diag(ones(1, cub-sz_x*sz_y), sz_x*sz_y)/dz ;
D_3d((1:sz_x*sz_y)*sz_y, (1:sz_x*sz_y)*sz_y + 1) = 0 ;
D_3d = D_3d(1:cub, 1:cub) ;
   

D2_3d = -2*diag(ones(1, cub), 0)/dx^2 - 2*diag(ones(1, cub), 0)/dy^2 -...
    2*diag(ones(1, cub), 0)/dz^2 + diag(ones(1, cub-1), 1)/dx^2 + ...
    diag(ones(1, cub-1), -1)/dx^2 + diag(ones(1, cub-sz_x), -sz_x)/dy^2 + ...
    diag(ones(1, cub-sz_x), sz_x)/dy^2 + diag(ones(1, cub-sz_x*sz_y), sz_x*sz_y)/dz^2 +...
    diag(ones(1, cub-sz_x*sz_y), -sz_x*sz_y)/dz^2 ;

D2_3d((1:sz_x*sz_y)*sz_y + 1, (1:sz_x*sz_y)*sz_y) = 0 ;
D2_3d((1:sz_x*sz_y)*sz_y, (1:sz_x*sz_y)*sz_y + 1) = 0 ;
D2_3d = D_3d(1:cub, 1:cub) ;

%% Integrating  
r2_min_r1 = sqrt((x_s-x_m).^2 + (y_s-y_m).^2 + (z_s-z_m).^2) ;
X = X(:) ; Y = Y(:) ; Z = Z(:) ;
r2_min_r = sqrt((X-x_m).^2 + (Y-y_m).^2 + (Z-z_m).^2) ;
r1_min_r = sqrt((X-x_s).^2 + (Y-y_s).^2 + (Z-z_s).^2) ;    

figure;
hold all ;


    %% Calculating Mach number
    Mxy = Mach(X, Y, Z, V, d, dz) ;
    Mxy_mat = repmat(Mxy.', cub, 1) ;
    P = f_array*0 ;
    %for each frequency 
    for f_index = 1:length(f_array)
        f = f_array(f_index) ;
        omega = 2*pi*f ;
        k = omega/c0 ;
        % Building the kernel of the equation
        Matr = -2i*Mxy_mat*omega.*D_3d/c0 + Mxy_mat.^2.*D2_3d ;
        p_i = 1/32/pi^4*exp(1i*k*r1_min_r)./r1_min_r./(Mxy.^2 - 1) ;
        p_ip1 = X ;
        for i=1:8
            if (i>1) p_i = p_ip1 ; end ;
            p_ip1 = 1/(4*pi)*sum(Matr * (p_i .* exp(1i*k*r2_min_r)./r2_min_r)...
                *dx *dy *dz ) ;       
        end
    end %end for each frequency

%     figure; 
%     plot(f_array, P) ;
    
    plot((1:length(f_array))*c0/f_array(end), real(fft(p_ip1))) ;

%     save(['lippman_data_', num2str(V)], 'P') ;

