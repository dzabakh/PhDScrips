function outvar = fb_contur(k_ar,h_s,h_m,L,rho, theta_exp, j)
    k = k_ar(end);
    dk = k/5000;
    
    
    k_re1 = [dk:dk:k] ;
    k_re2 = [(k+dk):dk:2*k] ;
    k_re = [k_re1, k_re2] ;
    
    k_im = exp((k_re2 - k)/10) - 1;
    k_cont1 = [k_re1, (k_re2 + 1i*k_im)] - 1i * gaussmf(k_re, [k/30 k]);
    k_cont2 = [k_re1, (k_re2 - 1i*k_im)] - 1i * gaussmf(k_re, [k/30 k]);
    
    dk_ar1 = abs(k_cont1(2:end) - k_cont1(1:end-1));   
    dk_ar2 = abs(k_cont2(2:end) - k_cont2(1:end-1));   
    
%     k_par = [0:dk:k-dk] ;
    
    theta2 = theta_exp ;
    theta2(j) = [] ;
    l_j1 = ones(size(k_re)) ; 
    l_j2 = ones(size(k_re)) ; 
    
%     dk_ar = k_par(2:end) - k_par(1:end-1);

    for i = 1:(length(theta_exp) - 1)
        l_j1 = l_j1 .* (asin(k_cont1/k) - theta2(i))/(theta_exp(j) - theta2(i)) ;
        l_j2 = l_j2 .* (asin(k_cont2/k) - theta2(i))/(theta_exp(j) - theta2(i)) ;
    end  

    P1 = (1i*rho/(8*pi)) .* exp(1i * (h_s + h_m) * my_sqrt((k^2 - k_cont1.^2),k)) .* besselh(0,1,k_cont1*L) .* k_cont1 ./ my_sqrt((k^2 - k_cont1.^2),k) .* l_j1;
    P2 = (1i*rho/(8*pi)) .* exp(1i * (h_s + h_m) * my_sqrt((k^2 - k_cont2.^2),k)) .* besselh(0,2,k_cont2*L) .* k_cont2 ./ my_sqrt((k^2 - k_cont2.^2),k) .* l_j2;
        
    P1(isnan(P1)) = 0 ;
    P2(isnan(P2)) = 0 ;
    P1 = (P1(1:end-1) + P1(2:end)) /2;
    P2 = (P2(1:end-1) + P2(2:end)) /2;
    P = P1*dk_ar1.' + P2*dk_ar2.' ;
outvar = P ;