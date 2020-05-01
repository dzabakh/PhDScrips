% load('U35mps') ;
% V_prof = U_m_40cm(7:end,1) ;
% 
% M = V_prof*V_prof' ;
% M = M(:) ;

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

dist_source_micro = sqrt((x_source-x_micro)^2 + (y_source-y_micro)^2) ; 
dist_source_stream = sqrt((X-x_source).^2 + (Y-y_source).^2) ; 
dist_micro_stream  = sqrt((X-x_micro ).^2 + (Y-y_micro ).^2) ; 

% R = sqrt(X.^2 + Y.^2) ;

% M = V/c0 * 0.5 * (1 - tanh((R - 0.2)/0.02)) ;
M = V/c0 * 0.5 .* (1 - tanh((dist - 0.290)/0.02)) ;


results_ar = [] ; 

f_ar = (200:200:5000) ;


for f =   f_ar
    
    f
    
    omega = 2*pi*f ; 
    k = omega / c0 ;
    
    result = 0 ; 
    
    res_ar = [];
    
    k_step =  0.02 ; 
    
    % beta = k_z
    for beta = 2*pi*(-1/step/4 : k_step : 1/step/4) ; 
        % omega^2/c^2 - k_z^2
        lambda = sqrt(k^2 - beta^2) ;
        
        % distances array
        r_array = (step : step/10 : 3) ; 
        % Green's function
        wave = -1i/4 * besselh(0,1,lambda*r_array) ;
        
        % adding zero to r_array and to Green's function
        r_array = [0, r_array] ;
        wave = [0 , wave] ;
        
        % interpolate Green's function to dist
        matrix = interp1(r_array , wave , dist) ;
        
        % dist_source_stream - distance to source in (x,y) plane
        % interpolate Green's function to dist_source_stream
        psi0 = interp1(r_array , wave , dist_source_stream) ;
        % yadro intura
        coeff = (- M.^2 * beta^2 + 2 * M * k * beta) * step^2 ; 
        % vsyo krome yadra intura
        reconst = interp1(r_array , wave , dist_micro_stream) .* coeff  ;
        
        mat_coef = - matrix .* (ones(sz,1) * coeff.')  + eye(sz,sz); 
        
        psi = mat_coef \ psi0 ;
        
        reply = (interp1(r_array , wave , dist_source_micro) + ...
            reconst.' * psi) * exp(1i*beta*(z_micro-z_source)) ; 
        
        result = result + reply / (2*pi) * k_step *2*pi; 
       
        res_ar = [res_ar , reply] ; 
        
        
    end
    
    result
    
    results_ar = [results_ar , result] ; 
        
end

figure ;
plot(results_ar);


% figure
% hold on
% 
% surf(X,Y,M) 
    


