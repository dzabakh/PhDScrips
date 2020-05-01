x_source = -0.3 ;
y_source = 0 ; 
z_source = 1 ; 

x_micro = 0.3 ; 
y_micro = 0 ; 
z_micro = 0.4 ; 

step = 0.02;

c0 = 343 ; 

V_array = [0, 20, 40, 60, 80] ;
for v_ind = 1:5
V = V_array(v_ind)

x_array = (-0.24:step:0.24) ;
y_array = (-0.24:step:0.24) ;

[X,Y]  = meshgrid(x_array, y_array) ;
X = X(:) ;
Y = Y(:) ;
sz = length(X) ; 

dist_x = X * ones(1,sz) - ones(sz,1)*X' ;
dist_y = Y * ones(1,sz) - ones(sz,1)*Y' ;
dist = sqrt(dist_x.^2 + dist_y.^2) ;

dist_source_micro = sqrt((x_source-x_micro)^2 + (y_source-y_micro)^2) ; 
dist_source_stream = sqrt((X-x_source).^2 + (Y-y_source).^2) ; 
dist_micro_stream  = sqrt((X-x_micro ).^2 + (Y-y_micro ).^2) ; 

R = sqrt(X.^2 + Y.^2) ;

M = V/c0 * 0.5 * (1 - tanh((R - 0.2)/0.02)) ;
% M = V/c0 * 0.5 .* (1 - tanh((Y - 0.290)/0.020)) ;


results_ar = [] ; 

f_ar = (200:200:3000) ;


for f =   f_ar
    
    f
    
    omega = 2*pi*f ; 
    k = omega / c0 ;
    
    result = 0 ; 
    
    res_ar = [];
    
    k_step =  0.02 ; 
    
     
    for beta = 2*pi*(-1/step/4 : k_step : 1/step/4) ; 
        lambda = sqrt(k^2 - beta^2) ;
     
        r_array = (step : step/10 : 3) ; 
        
        wave = -1i/4 * besselh(0,1,lambda*r_array) ;
        
        r_array = [0, r_array] ;
        wave = [0 , wave] ;
        matrix = interp1(r_array , wave , dist) ;
        
        psi0 = interp1(r_array , wave , dist_source_stream) ;
        
        coeff = (- M.^2 * beta^2 + 2 * M * k * beta) * step^2 ; 
        
        reconst = interp1(r_array , wave , dist_micro_stream) .* coeff  ;
        
        mat_coef = - matrix .* (ones(sz,1) * coeff.')  + eye(sz,sz); 
        
        psi = mat_coef \ psi0 ;
        
        reply = (interp1(r_array , wave , dist_source_micro) + reconst.' * psi) * exp(1i*beta*(z_micro-z_source)) ; 
        
        result = result + reply / (2*pi) * k_step *2*pi; 
       
        res_ar = [res_ar , reply] ; 
        
        
    end
    
    result
    
    results_ar = [results_ar , result] ; 
        
end

figure ;
plot(results_ar);
save(['results_lshV', num2str(V)], 'results_ar', 'f_ar') ;

end
