% Function Kij integrates the Kernel with one of two shape functions on 
% the interval [x1,x2]
% If flag == 0 the The kernel is integrated with ascending shape function 
%(x_ar - x1)/(x2-x1). In the case flag == 1 it is integrated with descending
% function  

function [res]= Kij_full(k, x0, x2, x1, phi_ar, phi0, alpha, n, n_im, imag, flag)

    % k is a wavenumber
    % a is tangent of the cone angle
    % n is the number of points of discretization on the real axis
    % 2*n_im is the is the number of points of discretization on the 
    % imaginary axis

    ka2 = k * alpha^2 ;
    dx = x2/n ;
    dx_im = imag/n_im ;


    %  Integration contours is introduced.  The 
    % contour is deformed into the upper half plane.

    x_ar_sing_1 = x1 + 1i*(dx_im:dx_im:imag) ;
    x_ar_sing_2 = 1i * imag + (x1 + dx : dx : x2) ;
    x_ar_sing_3 = x2 + 1i * (imag - dx_im : -dx_im : dx_im) ;

    x = [x_ar_sing_1, x_ar_sing_2, x_ar_sing_3] ;
    dx_ar_sing = x(2 : end) - x(1 : end-1) ;

    ones_phi = ones(length(phi_ar), 1) ;

    x_diff_grid = ones_phi .* (x0 - x) ;
    x_grid = ones_phi * x ;

    phi_grid = (phi_ar.' - phi0) * ones(1,length(x)) ;
    dphi_ar = phi_ar(2:end) - phi_ar(1:end-1);

    % Here one of the shape functions is chosen 

    if (flag==0)
        N_sing = (x_grid - x1) / (x2 - x1) ;
    elseif (flag==1)    
        N_sing = (x2 - x_grid) / (x2 - x1) ;
    end 

    % Kernel is calculated

    Kernel = 1i * x_grid .* ka2 / (2 * pi) .* (1 ./ x_diff_grid + (x_grid - x0 * cos(phi_grid)) ./ x_diff_grid .^ 2)...
            .* exp(0.5i * ka2 .* (x0^2 + x_grid.^2 - 2 * x_grid .* x0 .* cos(phi_grid)) ./ x_diff_grid) .* N_sing ;

    Kernel = (Kernel(2:end, :) + Kernel(1:end-1, :))/2 ;
    Kernel = (Kernel(:, 2:end) + Kernel(:, 1:end-1))/2 ;
    Kernel = Kernel .* ((dphi_ar.') * dx_ar_sing) ;


    %%% integral is evaluated 

    res = sum(Kernel,2);

end