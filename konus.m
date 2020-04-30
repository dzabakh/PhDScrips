clear all;

%% Loading experimental data
folder1_2 = '2018.5.10/cal10cm_big/'; % calibration experiment for mic 2 (KC)
folder1_4 = '2017.6.9/10cm_4/'; % calibration experiment for mic 4
folder1_6 = '2017.6.9/10cm_6/'; % calibration experiment for mic 6
folder2 = '2017.7.20/pricel/'; % positioning experiment
folder3 = '2018.9.8/ogurets_test/'; % experiment with cone

chan12_cal = wavread([folder1_2, 'chan1-2.wav']);
chan34_cal = wavread([folder1_4, 'chan3-4.wav']);
chan56_cal = wavread([folder1_6, 'chan5-6.wav']);

chan34_pos = wavread([folder2, 'chan3-4.wav']);
chan56_pos = wavread([folder2, 'chan5-6.wav']);

chan12_exp1 = wavread([folder3, 'chan1-2.wav']);

[mls, Fs] = wavread('mls_long');
mls = mls(:,1) ;
sz_mls = length(mls) ;

%% Syncing experimental data
xcorr12_cal = real(fft(ifft(chan12_cal(1:sz_mls,1)) .* conj(ifft(mls)))) ;
xcorr34_cal = real(fft(ifft(chan34_cal(1:sz_mls,1)) .* conj(ifft(mls)))) ;
xcorr56_cal = real(fft(ifft(chan56_cal(1:sz_mls,1)) .* conj(ifft(mls)))) ;

xcorr12_exp1 = real(fft(ifft(chan12_exp1(1:sz_mls,1)) .* conj(ifft(mls)))) ;
xcorr34_pos = real(fft(ifft(chan34_pos(1:sz_mls,1)) .* conj(ifft(mls)))) ;
xcorr56_pos = real(fft(ifft(chan56_pos(1:sz_mls,1)) .* conj(ifft(mls)))) ;

[max12_cal, t_start12_cal] = max(xcorr12_cal) ;
[max34_cal, t_start34_cal] = max(xcorr34_cal) ;
[max56_cal, t_start56_cal] = max(xcorr56_cal) ;
[max12_exp1, t_start12_exp1] = max(xcorr12_exp1) ;
[max34_pos, t_start34_pos] = max(xcorr34_pos) ;
[max56_pos, t_start56_pos] = max(xcorr56_pos) ;

t_end12_cal = t_start12_cal + sz_mls - 1 ;
t_end34_cal = t_start34_cal + sz_mls - 1 ;
t_end56_cal = t_start56_cal + sz_mls - 1 ;
t_end12_exp1 = t_start12_exp1 + sz_mls - 1 ;
t_end34_pos = t_start34_pos + sz_mls - 1 ;
t_end56_pos = t_start56_pos + sz_mls - 1 ;

pos1 = chan34_pos(t_start34_pos:t_end34_pos, 2).' / 31.6;
pos2 = chan56_pos(t_start56_pos:t_end56_pos, 2).' / 31.6;
cal_pos1 = chan34_cal(t_start34_cal:t_end34_cal, 2).' / 31.6;
cal_pos2 = chan56_cal(t_start56_cal:t_end56_cal, 2).' / 31.6;

adapter1 = chan12_cal(t_start12_cal:t_end12_cal, 2).' / 31.6;
probemic = chan12_exp1(t_start12_exp1:t_end12_exp1, 2).' / 31.6;

%% Correlation, MLS method

sp_pos1 = ifft(pos1) ;
sp_pos2 = ifft(pos2) ;
sp_cal_pos1 = ifft(cal_pos1) ;
sp_cal_pos2 = ifft(cal_pos2) ;
sp_a1 = ifft(adapter1) ; 
sp_probe1 = ifft(probemic) ;
sp_mls = ifft(mls.') ;

corr_pos1 = real(fft(sp_pos1 .* conj(sp_mls))) ;
corr_pos2 = real(fft(sp_pos2 .* conj(sp_mls))) ;
corr_cal_pos1 = real(fft(sp_cal_pos1 .* conj(sp_mls))) ;
corr_cal_pos2 = real(fft(sp_cal_pos2 .* conj(sp_mls))) ;
corr_a1 = real(fft(sp_a1 .* conj(sp_mls))) ;
corr_probe = real(fft(sp_probe1 .* conj(sp_mls))) ;

sz = 1000 ;
c = 343;
rho = 1.21 ;
f_max = Fs ; 
f_array = (0:sz-1)*f_max /sz ; 
df = f_max/sz ;
k_source = 2*pi*f_array/c;

corr_a1 = corr_a1(1:sz) ;
corr_probe = corr_probe(1:sz) ;
corr_pos1 = corr_pos1(1:sz) ;
corr_pos2 = corr_pos2(1:sz) ;
corr_cal_pos1 = corr_cal_pos1(1:sz) ;
corr_cal_pos2 = corr_cal_pos2(1:sz) ;

sp_corr_a1 = ifft(corr_a1);
sp_corr_probe = ifft(corr_probe);
sp_corr_pos1 = ifft(corr_pos1);
sp_corr_pos2 = ifft(corr_pos2);
sp_corr_cal_pos1 = ifft(corr_cal_pos1);
sp_corr_cal_pos2 = ifft(corr_cal_pos2);

f_high = 10000 ; 
f_width =  200 ; 
    
mask_1 = 0.5*(1-tanh((f_array - f_high)/f_width)) + 0.5*(1+tanh((-f_max + f_array + f_high)/f_width)) ; 

    
f_low = 500 ; 
f_width = 200 ; 

mask_2 = (0.5*(1+tanh((f_array - f_low)/f_width))) .*( 0.5*(1-tanh((-f_max + f_array + f_low)/f_width))) ; 

mask = mask_1.* mask_2 ; 
mask(length(mask)/2 + 1: end) = 0; 


%% Calibration 

delta = 0 * mask ; 
delta(1) = 1; 
delta_freq = fft(delta) ; 
delta  = real(ifft(mask.*delta_freq)) ;
normirovka = delta(1) ;
mask = mask/normirovka;

%% Positioning 

sp_corr_pos1 = (1/sz) * (4*pi/rho) * sp_corr_pos1./ sp_corr_cal_pos1.*exp(1i*k_source*0.1); 
sp_corr_pos1 = sp_corr_pos1 .* mask;
sig_pos1 = abs(fft(sp_corr_pos1)) ; 

sp_corr_pos2 = (1/sz) * (4*pi/rho) * sp_corr_pos2./ sp_corr_cal_pos2.*exp(1i*k_source*0.1); 
sp_corr_pos2 = sp_corr_pos2 .* mask;
sig_pos2 = abs(fft(sp_corr_pos2)) ; 

   
% figure;
% plot((1:sz)*c/Fs,sig_pos1);
% hold all; 
% plot((1:sz)*c/Fs,sig_pos2);



%% Building impulse response
sp_corr_probe = (1/sz) * (4*pi/rho) * sp_corr_probe./ sp_corr_a1.*exp(1i*k_source*0.1); 
sp_corr_probe = sp_corr_probe .* mask;
sig_probe = real(fft(sp_corr_probe)) ; 


figure ; 
plot((1:sz)*c/Fs, sig_probe);
hold all; 

%% Working with a piece of impulse response
 
r = 1.2 ;
r_low   = r - 0.4 ;
r_high  = r + 0.4 ;
r_width = 0.02 ;

t_low   = r_low   * Fs/c ;
t_high  = r_high  * Fs/c ;
t_width = r_width * Fs/c ;

mask_t1 = 0.5*(1-tanh(((0:sz-1) - t_high)/t_width));
mask_t2 = 0.5*(1+tanh(((0:sz-1) - t_low)/t_width)) ;
 
mask_t = mask_t1 .* mask_t2 ;
 
plot((1:sz)*c/Fs, mask_t) ;

sig_interest = sig_probe.* mask_t ;

sp_sig_interest = ifft(sig_interest) ;
figure; plot(f_array, real(sp_sig_interest)) ;
