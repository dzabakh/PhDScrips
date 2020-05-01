 load data_01 % свободное поле
% load data_07 % с экраном, микрофон на оси
% load data_08 % с экраном, микрофон на 20 см в освещенной зоне
% load data_09 % с экраном, микрофон на 40 см в освещенной зоне
% load data_10 % с экраном, микрофон на 20 см в тени
%load data_11 % с экраном, микрофон на 40 см в тени






[MLS_sample, Fs] = wavread('MLS');



dt_rec = (chan1x(end)-chan1x(1)) /  (length(chan1x)-1) ;
dt_sam = 1/Fs ; 

sz_sam = length(MLS_sample);


t_array = (1:500000)*1e-5 ; 
nn = floor(1e5*dt_sam * sz_sam)-3 ; 
t_array_1 = t_array(1:nn) ;

%%%%%%%%%%%%%%%% interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig_1 =interp1(chan1x , chan1y , t_array) ; 
sig_2 =interp1(chan2x , chan2y , t_array) ; 
mls_1   =interp1((0:sz_sam-1)*dt_sam, MLS_sample , t_array_1) ;  
mls = 0*sig_1 ; 
mls(1:nn) = mls_1 ; 

%%%%%%%%%%%%%%%% correlation processing %%%%%%%%%%%%%%%%%%%%%%%%%

sp_1 = ifft(sig_1) ; 
sp_2 = ifft(sig_2) ;
sp_mls = ifft(mls) ; 

corr_1 = real(fft(sp_1 .* conj(sp_mls))) ;
corr_2 = real(fft(sp_2 .* conj(sp_mls))) ;

%%%%%%%%%%%%%%%% rought time mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dum, idx_max] = max(corr_1) ; 

idx_max = idx_max - 200 ; 
start_point = mod(idx_max , 500000) ; 

if start_point > 0 
    corr_1 = [ corr_1(start_point+1 : end) , corr_1(1:start_point) ] ;
    corr_2 = [ corr_2(start_point+1 : end) , corr_2(1:start_point) ] ;
end

corr_1 = corr_1(1:3000) ; 
corr_2 = corr_2(1:3000) ;

% idx_active = (410501:413500) ; %file 01
%idx_active = (497001:500000) ; %file 09
%corr_1 = corr_1(idx_active) ; 
%corr_2 = corr_2(idx_active) ;

% figure
% plot(corr_1)
% hold all
% plot(corr_2)

sz = 3000 ; 


%%%%%%%%%%%%%%%%% compressing by the simplest formula %%%%%%%%%%%%%%

c = 330;
f_max = 1e5 ; 
f_array = (0:sz-1)*f_max / sz ; 
R_tube = 0.01875; % radius of the tube
R = 0.05 ; 
rho = 1.29 ; 
L_ref = 1 ; 

sz_short = floor(12000 / (f_array(2)-f_array(1))) ; 
k_source = 2*pi*f_array(1:sz_short)/c;
sp_corr_1 = ifft(corr_1) ; 
sp_corr_1_s = sp_corr_1(1:sz_short) ; 
% piston source
% W_w_p = 2*k_source.*(pi*R_tube^2/(2*rho)).*sp_corr_1_s.*exp(-(1i*k_source/2)*(sqrt(R^2 + R_tube^2) + R))...
%     ./sin((k_source/2)*(sqrt(R^2 + R_tube^2) - R));
% W_w_p(1) = W_w_p(2);
W_w_m = exp(-1i*k_source*R).*sp_corr_1_s*(4*pi*R/rho);
%save calibrator f_array sp_corr_1_s
%load calibrator f_array sp_corr_1_s


sp_corr_2 = ifft(corr_2) ; 
sp_corr_2_s = sp_corr_2(1:sz_short) ; 

  sp_corr_2_red = (4*pi/rho)*sp_corr_2_s ./ W_w_m; %correction using monopole source 
  %sp_corr_2_red = sp_corr_2_s ./W_w_p ; %correction using piston source
sp_corr_2 = 0*sp_corr_2 ; 

sp_corr_2(1:sz_short) = sp_corr_2_red ; 

sig_2_red = 2*real(fft(sp_corr_2))/sz ; 




%%%%%%%%%%%%%%%% filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_high = 4000 ; 
f_width = 1500 ; 

mask_1 = 0.5*(1-tanh((f_array - f_high)/f_width)) + 0.5*(1+tanh((-f_max + f_array + f_high)/f_width)) ; 

f_low = 200 ; 
f_width = 70 ; 

mask_2 = (0.5*(1+tanh((f_array - f_low)/f_width))) .*( 0.5*(1-tanh((-f_max + f_array + f_low)/f_width))) ; 

mask = mask_1 .* mask_2 ; 

    testarray    =  zeros(1,sz);
    testarray(1) = 1; 
    testarray_freq = fft(testarray) ; 
    testarray  = real(ifft(mask.*testarray_freq)) ;
    maxval = testarray(1) ;
mask = mask/maxval;
%sp_corr_1 = ifft(corr_1) ; 
sp_corr_2 = ifft(sig_2_red) ; 

%corr_1 = real(fft(sp_corr_1 .* mask)) ; 
sig_2_red = real(fft(sp_corr_2 .* mask)) ; 

figure
plot((1:sz)*1e-5*343 , sig_2_red)





