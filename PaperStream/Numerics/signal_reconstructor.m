figure;
hold all ;
V = 0 ;
load(['results_lshV', num2str(V),'.mat']) ;


f_ar_new = (0:200:39900) ;

sp = 0* f_ar_new ;

%results_ar = -exp(1i*2*pi*f_ar / 343 *1) / (4 * pi * 1) ;

recon_factor = exp(-1i*2*pi*f_ar / 343 *0.6) ;


sp(2:16) = results_ar .* recon_factor ;

f_high = 3000;
    f_width = 200;

    mask_1 = 0.5*(1-tanh((f_ar_new - f_high)/f_width)) + 0.5*(1+tanh((-f_max + f_ar_new + f_high)/f_width)) ; 

    f_low = 200 ; 
    f_width = 200 ; 

    mask_2 = (0.5*(1+tanh((f_ar_new - f_low)/f_width))) .*( 0.5*(1-tanh((-f_max + f_ar_new + f_low)/f_width))) ; 

    mask = mask_1 .* mask_2 ;


sp = sp.* mask ; 



signal = 2*real(fft(sp)) * 200 / (2*pi) /85;

t_array = (1:200)*0.25e-4 ;


plot(t_array*343 + 0.6161 , -2*signal/1.95*1.284)

