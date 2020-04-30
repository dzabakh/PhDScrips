clear all;

 %  folder = '2019.6.13\empty_room_big_speaker\';
 folder = '2019.6.13\2_off_room_big_speaker\';
chan12 = audioread([folder, 'chan1-2.wav'],[1,11*48000]);
chan34 = audioread([folder, 'chan3-4.wav'],[1,11*48000]);
chan56 = audioread([folder, 'chan5-6.wav'],[1,11*48000]);
chan78 = audioread([folder, 'chan7-8.wav'],[1,11*48000]);

load('cal_curves') ;

[mls, Fs] = audioread('mls_long.wav');
mls = mls(:,1) ;
sz_mls = length(mls) ;

xcorr12 = real(fft(ifft(chan12(1:sz_mls,1)) .* conj(ifft(mls)))) ;
xcorr34 = real(fft(ifft(chan34(1:sz_mls,1)) .* conj(ifft(mls)))) ;
xcorr56 = real(fft(ifft(chan56(1:sz_mls,1)) .* conj(ifft(mls)))) ;
xcorr78 = real(fft(ifft(chan78(1:sz_mls,1)) .* conj(ifft(mls)))) ;

[max12, t_start12] = max(xcorr12) ;
[max34, t_start34] = max(xcorr34) ;
[max56, t_start56] = max(xcorr56) ;
[max78, t_start78] = max(xcorr78) ;

t_end12 = t_start12 + sz_mls - 1 ;
t_end34 = t_start34 + sz_mls - 1 ;
t_end56 = t_start56 + sz_mls - 1 ;
t_end78 = t_start78 + sz_mls - 1 ;

adapter2 = chan56(t_start34:t_end34, 2).';
probemic = chan12(t_start12:t_end12, 2).';
bkmic = chan78(t_start78:t_end78, 2).';
adapter1 = chan34(t_start56:t_end56, 2).';%/ 31.6; %+30dB

sp_a1 = ifft(adapter1) ; 
sp_a2 = ifft(adapter2) ;
sp_probe = ifft(probemic) ;
sp_bkmic = ifft(bkmic) ;
sp_mls = ifft(mls.') ;
% sp_mls = ifft( chan34(t_start34:t_end34, 1).');
corr_a1 = real(fft(sp_a2 .* conj(sp_mls))) ;
corr_a2= real(fft(sp_a1 .* conj(sp_mls))) ;
corr_bkmic= real(fft(sp_bkmic .* conj(sp_mls))) ;
corr_probe = real(fft(sp_probe .* conj(sp_mls))) ;
 figure;
plot(corr_a1)
hold all
plot(corr_a2)
hold all
 plot(corr_probe)

sz = 3000 ;
c = 343;
R_tube = 0.035/2; % radius of the tube 
rho = 1.21 ; 
b = 0.09 ;
a = 0.07 ;
f_max = Fs ; 
f_array = (0:sz-1)*f_max /sz ; 
df = f_max/sz ;
k_source = 2*pi*f_array/c;

corr_a1 = corr_a1(1:sz) ;
corr_a2 = corr_a2(1:sz) ;
corr_bkmic = corr_bkmic(1:sz) ;
corr_probe = corr_probe(1:sz) ;




% chan12_cal = interp1(f_ar_cal,chan12_cal,f_array);
% chan56_cal = interp1(f_ar_cal,chan56_cal,f_array);
% chan12_cal(isnan(chan12_cal))=1;
% chan56_cal(isnan(chan56_cal))=1;
sp_corr_a2 = ifft(corr_a2);
sp_corr_a1 = ifft(corr_a1)./cal_a1;
sp_corr_probe = ifft(corr_probe)./cal_probe;
    

    
   
    




f_high = 4500 ; 
f_width =  500 ; 
    
mask_1 = 0.5*(1-tanh((f_array - f_high)/f_width)) + 0.5*(1+tanh((-f_max + f_array + f_high)/f_width)) ; 

    
f_low = 200 ; 
f_width = 200 ; 

mask_2 = (0.5*(1+tanh((f_array - f_low)/f_width))) .*( 0.5*(1-tanh((-f_max + f_array + f_low)/f_width))) ; 

mask = mask_1;%.* mask_2 ; 
mask(length(mask)/2 + 1: end) = 0; 


  %%%%%%%%%%%%%%%%%%%%%Calibration%%%%%%%%%%%%%%%%%%%

    delta = 0 * mask ; 
    delta(1) = 1; 
    delta_freq = fft(delta) ; 
    delta  = real(ifft(mask.*delta_freq)) ;
    normirovka = delta(1) ;
    mask = mask/normirovka;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

sp_corr_probe = sp_corr_probe .* mask;
   
    

%A = (-sp_corr_a2 .* exp(1i*k_source*b) + sp_corr_a1 .* exp(1i*k_source*a)) ./ (exp(2i*k_source*a) - exp(2i*k_source*b)) ;
A = (sp_corr_a2 .* exp(1i*k_source*(2*a+b)) - sp_corr_a1 .* exp(1i*k_source*(a+2*b))) ./ (exp(2i*k_source*a) - exp(2i*k_source*b)) ;  
W = k_source.*A*2i*pi*R_tube^2/(rho);
W(1) = W(2) ;     
    
sp_corr_probe = (1/sz) * (4*pi/rho) * sp_corr_probe ./ W; 

    
sig_probe = -real(fft(sp_corr_probe)) ; 
[amp,idx] = max(sig_probe);
idx = idx*c/Fs;
% norm = idx*amp ;
figure ; 
% sig_probe = sig_probe; %/norm;
plot((1:sz)*c/Fs, sig_probe);
xlabel('$ct$, m','Interpreter','latex','FontSize',14)
ylabel('Амплитуда ','Interpreter','latex','FontSize',14)
% hold all;
% % % 
% % % 
% % % %%%%%Working with a piece of impulse response%%%%%%%%%%%%
% % % 
% r = 1.07 ;
% r_low   = r - 0.17 ;
% r_high  = r + 0.17 ;
% r_width = 0.02 ;
% 
% t_low   = r_low   * Fs/c ;
% t_high  = r_high  * Fs/c ;
% t_width = r_width * Fs/c ;
% 
% mask_t1 = 0.5*(1-tanh(((0:sz-1) - t_high)/t_width));
% mask_t2 = 0.5*(1+tanh(((0:sz-1) - t_low)/t_width)) ;
%  
% mask_t = mask_t1 .* mask_t2 ;
%  
% plot((1:sz)*c/Fs, mask_t) ;
% 
% sig_interest = sig_probe.* mask_t ;
% figure;
% plot((1:sz)*c/Fs,sig_interest) ;
%   
% sig_probe_fs = sig_probe ;
% R_fs0 = ifft(sig_interest) ;
% nfft = length(R_fs0) ;
% R_fs0(nfft/2+1:end) = 0 ;
% R_fs0 = R_fs0*2 ;

% save('R_fs0', 'R_fs0', 'sz' , 'Fs', 'sig_probe_fs', 'f_array') ;

% sig_probe_wall = sig_probe ;
% R_wall0 = ifft(sig_interest) ;
% nfft = length(R_wall0) ;
% R_wall0(nfft/2+1:end) = 0 ;
% R_wall0 = R_wall0*2 ;
% 
% save('R_wall0', 'R_wall0', 'sz' , 'Fs', 'sig_probe_wall', 'f_array') ;

