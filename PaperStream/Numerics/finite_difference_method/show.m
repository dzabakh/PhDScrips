clear all ;

V_array = [0,20,40,60,80] ;
figure; hold all;
for i = 1:5

V = V_array(i) ;
load(['res3d50_', num2str(V), '.mat'])

c = 343 ;
t = [dt:dt:n_t*dt] ;

Fs = 1/dt ;

f_array = (0:n_t-1)*Fs /n_t ; 
% f_high = 3000 ; 
% f_width = 200 ; 
f_high = 1500 ; 
f_width = 200 ; 
f_max = Fs ; 
mask_1 = 0.5*(1-tanh((f_array - f_high)/f_width)) + ...
         0.5*(1+tanh((-f_max + f_array + f_high)/f_width)) ; 

f_low = 100 ; 
f_width = 70 ; 
% f_low = 100 ; 
% f_width = 200 ; 
mask_2 = (0.5*(1+tanh((f_array - f_low)/f_width))) ...
       .*( 0.5*(1-tanh((-f_max + f_array + f_low)/f_width))) ; 
mask = mask_1 .* mask_2 ; 

delta = 0* mask ; 
delta(1) = 1 ; 
delta_filtered = real(fft(ifft(delta) .* mask)); 
normirovka = delta_filtered(1) ; 
mask = mask/normirovka ;

sp = mask.*ifft(res_mic) ; 

corr_f = real(fft(sp)) ;
% plot(1.655+c*dt*(1:n_t),-real(corr_f(1:n_t))/7.422e-6*0.3) ;
plot(c*dt*(1:n_t),corr_f) ;
end
