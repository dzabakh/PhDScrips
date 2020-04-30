clear all
[MLS_sample, Fs] = wavread('D:/ImpedanceTube/Processing/mls_long.wav');
MLS_sample = MLS_sample(:,2) ;
sz_sam = length(MLS_sample);
    

[mic_sig, Fs_rec] = wavread('D:/ImpedanceTube/Processing/2018.3.15/open1/chan3-4.wav') ;

sz_sam = length(MLS_sample);

%% Choosing the start moment
xcorr_mic = fft(ifft(mic_sig(1:sz_sam,2)) .* conj(ifft(MLS_sample))) ;
[maxval, t_start] = max(abs(xcorr_mic)) ;
t_start = t_start - 10 ;
t_end = t_start + sz_sam - 1 ;
mic_sig = mic_sig(t_start:t_end, 2).' ;

%% Corellation
sp_mic = ifft(mic_sig) ;
sp_mls = ifft(MLS_sample) ;
corr_mic = real(fft(sp_mic.' .* conj(sp_mls))) ;


sz = 1000 ;
c = 343;
rho = 1.21 ;
f_max = Fs ; 
f_array = (0:sz-1)*f_max /sz ; 
df = f_max/sz ;
k_source = 2*pi*f_array/c;

sp_mic = sp_mic(1:sz) ;
corr_mic = corr_mic(1:sz) ;
sp_corr_mic = ifft(corr_mic);


