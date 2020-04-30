clear all

c0 = 340 ;

L = 2 ;
H = 0.58 ;
h = 0.57 ;
r  = sqrt(L^2 + (H+h)^2) ;
rho_0 = 1.21 ;
%%%PART1%%%%%%%
ind = 0 ;
f = 50:50:6000 ;
for freq = f
    ind = ind + 1 ;
    k0 = freq*2*pi/340 ;
    omega = k0 * c0 ;
    dxi = 0.1 ;
    xi_start = floor((k0+10)/dxi)*dxi  ;
    xi_end   = floor((k0+100)/dxi)*dxi ;
    re_xi = xi_start : dxi : xi_end ;

    im_xi_1 = 1i*exp((re_xi-k0)/10) ;
    im_xi_1 = im_xi_1 - im_xi_1(1) ;
    xi_1 = re_xi + im_xi_1 ;
    xi_1 = [(0:dxi:k0+10),xi_1]  ;

    g = 2 * gaussmf(real(xi_1), [k0/5 k0]);
    xi_1 = xi_1 - 1i*g ; 
%      figure ; 
%      plot(real(xi_1), imag(xi_1));


    koren = my_sqrt((k0^2 - xi_1.^2), k0) ;
    integrand_1 = xi_1 .* besselh(0,1,xi_1 * L) ./ koren .* exp(1i * koren * (H + h)) .*biot_f(omega/2/pi, acos(xi_1/k0)) * dxi;
    int_1 = sum(integrand_1(2:end)) ;
%     figure; plot(abs(integrand_1)) ;


    %%%PART2%%%%%%%

    im_xi_2 = - 1i*exp((re_xi-k0)/10) ;
    im_xi_2 = im_xi_2 - im_xi_2(1) ;
    xi_2 = re_xi + im_xi_2 ;
    xi_2 = [(0:dxi:k0+10),xi_2] ;
    g = 2 * gaussmf(real(xi_2), [k0/5 k0]);
    xi_2 = xi_2 - 1i*g ; 
%      figure ; 
%      plot(real(xi_2), imag(xi_2));


    koren = my_sqrt((k0^2 - xi_2.^2), k0) ;
    integrand_2 = xi_2 .* besselh(0,2,xi_2 * L) ./ koren .* exp(1i * koren * (H + h)) .* biot_f(omega/2/pi, acos(xi_2/k0)) * dxi;
    int_2 = sum(integrand_2(2:end)) ;

%     figure; plot(abs(integrand_2));

     intsum(ind) = (int_1 + int_2) ;
% 
% 
% 
    end


f_max = 48000 ;
f_high = 4500 ; 
f_width =  500 ; 
    
mask_1 = 0.5*(1-tanh((f - f_high)/f_width)) + 0.5*(1+tanh((-f_max + f + f_high)/f_width)) ; 

    
f_low = 200 ; 
f_width = 100 ; 

mask_2 = (0.5*(1+tanh((f - f_low)/f_width))) .*( 0.5*(1-tanh((-f_max + f + f_low)/f_width))) ; 

mask = mask_1.* mask_2 ; 


fbd =  intsum * 1.*mask  ;
fb = rho_0 /2/pi * real(fft(2*pi*f.*fbd)) ;
figure;
plot(f,abs(fbd)) ;
