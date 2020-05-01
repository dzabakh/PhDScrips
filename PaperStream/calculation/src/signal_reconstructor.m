load data_calc_0

f_ar_new = (0:100:40000) ;


results_ar = p_sum ;
sp = 0* f_ar_new ;

% results_ar = -exp(1i*2*pi*f_array / 343 ) / (4 * pi * 1) ;

recon_factor = exp(-1i*2*pi*f_array / 343  * 0.965);


sp(2:31) = results_ar .* recon_factor ;



mask = 0.5* (1 - tanh((f_ar_new-4000)/1500)) ;
mask1= 0.5* (1 + tanh((f_ar_new-500)/150)) ;

mask = mask .* mask1  ; 

sp = sp.* mask ; 



signal = 2*real(fft(sp)) * 400 / (2*pi) /85;

t_array = (1:floor(length(f_ar_new)))*0.25e-4 ;

figure
plot(t_array*340 + 1-0.04 , signal , 'k')

