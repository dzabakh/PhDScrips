x_source = 1.62 ;
y_source = 0 ; 
z_source = 0 ; 

x_micro = -0.38 ; 
y_micro = 0.2 ; 
z_micro = -2.6 ; 

step = 0.02;

c0 = 343 ; 

V = 60 ; 

x_array = (-0.24:step:0.24) ;
y_array = (-0.24:step:0.24) ;

[X,Y]  = meshgrid(x_array, y_array) ;
X = X(:) ;
Y = Y(:) ;
sz = length(X) ; 

dist_x = X * ones(1,sz) - ones(sz,1)*X' ;
dist_y = Y * ones(1,sz) - ones(sz,1)*Y' ;
dist = sqrt(dist_x.^2 + dist_y.^2) ;

r2_min_r1 = sqrt((x_source-x_micro)^2 + (y_source-y_micro)^2) ; 
r_min_r2  = sqrt((X-x_micro ).^2 + (Y-y_micro ).^2) ;

M = V/c0 * 0.5 .* (1 - tanh((Y - 0.20)/0.02)) ;

results_ar = [] ; 

f_array = (200:200:5000) ;
res =[] ;

for f = f_array
    omega = 2*pi*f ; 
    k = omega / c0 ;
    result = 0 ; 
    k_step =  0.02 ; 

    for k_z = 2*pi*(-1/step/4 : k_step : 1/step/4) ;
        f
        k_z
        nev = 1000 ;
        k_xy = sqrt(k^2 - k_z^2) ;
        r_array = (step : step/10 : 3) ;
        wave = -1i/4 * besselh(0, 1, k_xy *r_array) ;
        r_array = [0, r_array] ;
        wave = [0, wave] ;
        G_xy_mat = interp1(r_array, wave, dist) ;    
        G_r2r1 = interp1(r_array, wave, r2_min_r1) ;
        G_r2r = interp1(r_array, wave, r_min_r2) ;
        kernel = (- M.^2 * k_z^2 + 2 * M * k * k_z) * step^2 ; 
        K = G_r2r .* kernel  ;

        f_x = exp(-1i *k_z *(z_micro - z_source)) * G_r2r1 ; 
        y_nmin1 = f_x * ones(sz,1) ;
        y_n = y_nmin1 * 0 ;
        
        while (abs(nev) > 1e-2)
           y_n = f_x + K * ones(1, sz) * y_nmin1 ;
           nev = sum(sqrt(y_nmin1.^2 - y_n.^2))/...
               sum(sqrt(y_nmin1.^2 + y_n.^2)) ;
           y_nmin1 = y_n ;
           y_n_num = sum(y_n) ;
        end
        result = result + y_n_num / (2*pi) * k_step *2*pi; 
    end    
    results_ar = [results_ar , result] ; 
end

figure;
plot(f_array, results_ar) ;

save('results_60') ;

