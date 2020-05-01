function [mask] = build_mask(f_array, f_high, f_low, f_low_width, f_high_width)

    f_max = f_array(end) ;
    mask_1 = 0.5*(1-tanh((f_array - f_high)/f_low_width)) + 0.5*(1+tanh((-f_max + f_array + f_high)/f_low_width)) ; 



    mask_2 = (0.5*(1+tanh((f_array - f_low)/f_high_width))) .*( 0.5*(1-tanh((-f_max + f_array + f_low)/f_high_width))) ; 

    mask = mask_1.* mask_2 ;
    mask(floor(end/2 + 1):end) = 0 ;
    

    delta = 0 * mask ; 
    delta(1) = 1; 
    delta_freq = fft(delta) ; 
    delta  = real(ifft(mask.*delta_freq)) ;
    normirovka = delta(1) ;
    mask = mask/normirovka;  

end    
    