clear all;
load('f_array') ;
angles_array = [0, 10, 20, 30, 40, 45, 50, 60, 70] * pi/ 180 ; % ot normali

L = 2 * sin(angles_array); %pi*(1/2-1/12)]) ;
H = 2 * cos(angles_array);%pi*(1/2-1/12)]) ;


for freq = 1:400
    f = f_array(freq) ; 

    om = 2*pi*f ; 
    c0 = 343 ; 
    k0 = om / c0 ; 
    rho0 = 1.2 ; 

    sin_angles_array = sin(angles_array) ; 
    sin_angles_array = [sin_angles_array , 1] ;

    ds = 0.001 ; 

    sin_array = [ds:ds:(1-ds)] ; 

    num_points = length(sin_angles_array) - 1 ; 

    shape_functions = zeros(length(sin_array) , num_points) ; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n = 1: num_points
       Z_desc = (sin_angles_array(n+1) - sin_array.') / (sin_angles_array(n+1) - sin_angles_array(n)) ; 
       Z_grow = (sin_array.' - sin_angles_array(n)) / (sin_angles_array(n+1) - sin_angles_array(n)) ; 

       idx = (sin_array >= sin_angles_array(n)) & (sin_array <= sin_angles_array(n+1)) ;

       shape_functions(idx , n)   = Z_desc(idx) ; 

       if n < (num_points - 0.5)
            shape_functions(idx , n+1) = Z_grow(idx) ; 
       end
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    num_experiments = length(L) ;


    bessel_matrix = besselj(0 , k0*(L.' * sin_array)) ; 

    gamma_matrix = k0 * ones(num_experiments , 1) * sqrt(1 - sin_array.^2) ; 
    
    exp_matrix = exp( 1i* k0 * ( H.' * sqrt(1 - sin_array.^2) )) ; 

    k_matrix = k0 * ones(num_experiments, 1) * sin_array ; 

    dk_matrix = ds * ones(1 , length(sin_array)) ; 
    dk_matrix(1) = dk_matrix(1) / 2 ; 
    dk_matrix(end) = dk_matrix(end) / 2 ; 
    dk_matrix = k0 * ones(num_experiments , 1) * dk_matrix ; 

    factor_matrix = (bessel_matrix ./ gamma_matrix) .* exp_matrix .* dk_matrix .* k_matrix; 

    M(:,:,freq) =  factor_matrix * shape_functions;
    M(isnan(M)) = 0;
    
    R = 2;%r1 + r2 ; 

    corrections = 1 ./( (-1i * om * rho0./ (4*pi*R)) .* exp(1i*k0*R) ) ;

    corrections = diag(corrections) ;
        
    matrix_corr(:,:,freq) = abs(corrections * M(:,:,freq) );
        
end
% save('M_matrix', 'M') ;


