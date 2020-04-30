clear
filename2 = 'D:\ImpedanceTube\Processing\2018.4.9\calibr\chan1-2.wav' ;
filename4 = 'D:\ImpedanceTube\Processing\2018.4.9\calibr\chan3-4.wav' ;
[sig2, Fs] = wavread(filename2) ;
sig4 = wavread(filename4) ;

mls_filename = 'D:\ImpedanceTube\Processing\2018.3.27\mls_long.wav' ;
[mls, Fs] = wavread(mls_filename) ;
mls = mls(:,2);
mls_len = length(mls) ;

ring2 = sig2(:,1) ;
ring4 = sig4(:,1) ;

corf2 = abs(fft(ifft(ring2(1:mls_len)).*conj(ifft(mls)))) ;
corf4 = abs(fft(ifft(ring4(1:mls_len)).*conj(ifft(mls)))) ;

[M, t_start2] = max(corf2) ;
[M, t_start4] = max(corf4) ;

t_end2 = t_start2 + mls_len -1;
mic2 = sig2(t_start2:t_end2, 2) ;

t_end4 = t_start4 + mls_len -1;
mic4 = sig4(t_start4:t_end4, 2) ;

dt = 1/Fs;
t_array = dt*(0:mls_len-1);
ct_array = 343*t_array;

corr_cal2 = (fft(ifft(mic2).*conj(ifft(mls))))';
corr_cal4 = (fft(ifft(mic4).*conj(ifft(mls))))';

mask = -tanh(5*(ct_array - 3));
mask = [mask(1:420),zeros(1,length(ct_array)-420)] ;
corr_cal2_cut = corr_cal2 .* mask ;
corr_cal4_cut = corr_cal4 .* mask ;

corr_cal2_cut = corr_cal2_cut(1:1000) ; 
corr_cal4_cut = corr_cal4_cut(1:1000) ; 

sz = length(corr_cal2_cut) ; 
f_array = (0:(sz-1)) / sz * Fs ; 

sp_2cal = ifft(corr_cal2_cut) ;
sp_4cal = ifft(corr_cal4_cut) ;

save('cal_curve', 'f_array', 'sp_2cal', 'sp_4cal') ;