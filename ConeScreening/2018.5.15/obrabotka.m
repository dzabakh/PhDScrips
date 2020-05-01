%% pricel
% 
% clear
% a = wavread('pricel/chan3-4.wav') ;
% a1 = a(:,1) ;
% a2 = a(:,2) ;
% b = wavread('pricel/chan5-6.wav') ;
% b1 = b(:,1) ;
% b2 = b(:,2) ;
% mls = wavread('mls_long') ;
% mls = mls(:,2) ;
% figure; 
% plot(real(fft(ifft(a1(1:length(mls))) .* conj(ifft(mls)))));
% hold all;  
% plot(real(fft(ifft(b1(1:length(mls))) .* conj(ifft(mls)))));
% a2=a2(1443:1442 + length(mls));
% b2=b2(1443:1442 + length(mls));
% figure; 
% plot(real(fft(ifft(a2) .* conj(ifft(mls)))));
% hold all; 
% plot(real(fft(ifft(b2) .* conj(ifft(mls)))));


%% obrabotka
clear
r = [1.858, 1.836, 1.822, 1.808, 1.808, 1.801, 1.815, 1.815, 1.836, 1.858, 1.879, 1.922, 1.951, 1.987, 2.037, 2.094] ;
mls = wavread('mls_long') ;
mls = mls(:,1) ;
N = length(mls) ;
[b, Fs] = wavread('tenj10/chan7-8.wav') ;

[m m_idx] = max(correl(b(:,1), mls)) ;

corr2 = correl(b(m_idx + 1: m_idx + N,2), mls ) ;

ct = (1:length(corr2)) * 343/ Fs ;

idx = (ct < 3) & (ct > 1);
corr2 = corr2(idx) ;
% hold all; 
% plot(ct(idx), corr2);
df = Fs/length(corr2) ;
f_array = df:df:Fs ;
hold all; plot(f_array, abs(ifft(corr2))* r(16))  ;
