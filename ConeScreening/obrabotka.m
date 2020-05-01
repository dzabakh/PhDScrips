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



clear

r = [1.858, 1.836, 1.822, 1.808, 1.808, 1.801, 1.815, 1.815, 1.836, 1.858, 1.879, 1.922, 1.951, 1.987, 2.037, 2.094] ;
[mls Fs] = wavread('mls_long') ;
mls = mls(:,1) ;
N = length(mls) ;    


%% kalibrovka
    cali_file = '2018.5.15/cal50/chan7-8.wav' ;
    cal = wavread(cali_file) ;
    [m m_idx] = max(correl(cal(:,1), mls)) ;
    corr_cal = correl(cal(m_idx + 1: m_idx + N, 2), mls) ;
    df = Fs/length(corr_cal) ;
    f_array = df:df:Fs ;
    k = 2*pi*f_array/343 ;
    
    sp_cal = ifft(corr_cal) .* exp(-1i*k(:)*0.45) * 0.45;
    ct = (1:length(corr_cal)) * 343/ Fs ;
    
    idx_cal = find((ct < 1.9)) ;
    ct = ct(idx_cal) ;
    corr_cal = fft(sp_cal) ;
    corr_cal = corr_cal(idx_cal) ;
    sp_cal = ifft(corr_cal) ;
    df_cal = Fs/length(sp_cal) ;
    f_cal = df_cal : df_cal :Fs ;
    
%% obrabotka
    [b, Fs] = wavread('2018.5.15/tenj10/chan7-8.wav') ;

    [m m_idx] = max(correl(b(:,1), mls)) ;

    corr2 = correl(b(m_idx + 1: m_idx + N, 2), mls ) ;
    ct = (1:length(corr2)) * 343/ Fs ;
%     figure; plot(ct, corr2) ;
    idx = (ct < 2.95) & (ct > 1.05);   
    
    corr2 = corr2(idx) ;
    corr2 = corr2(1:end-1) ;
    sp_corr2 = ifft(corr2) ./ sp_cal * r(5);

    corr2 = fft(sp_corr2) ;
    ct = (1:length(corr2)) * 343/ Fs + 1.05;
    
%     figure; plot(ct, real(corr2)) ;
    
    df = Fs/length(corr2) ;
    f_array = df:df:Fs ;
    mask = build_mask(f_array, 10000, 1000, 1000, 200) ;
    sp2 = ifft(corr2) .* mask(:) ;
%     figure; plot(f_array, sp2); hold all; plot(f_array, mask) ;
    
    hold all; plot(ct, real(fft(sp2))/length(sp2) ) ;
