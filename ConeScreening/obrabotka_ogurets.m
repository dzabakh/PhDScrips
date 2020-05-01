clear

r = [1.858, 1.836, 1.822, 1.808, 1.808, 1.801, 1.815, 1.815, 1.836, 1.858, 1.879, 1.922, 1.951, 1.987, 2.037, 2.094] ;
[mls Fs] = wavread('mls_long') ;
mls = mls(:,1) ;
N = length(mls) ;    


%% kalibrovka
    cali_file = '2018.9.8/cal40cm/chan1-2.wav' ;
    cal = wavread(cali_file) ;
    [m m_idx] = max(correl(cal(:,1), mls)) ;
    corr_cal = correl(cal(m_idx + 1: m_idx + N, 2), mls) ;
    df = Fs/length(corr_cal) ;
    f_array = df:df:Fs ;
    k = 2*pi*f_array/343 ;
    
    sp_cal = ifft(corr_cal) .* exp(-1i*k(:)*0.35) * 0.35;
    ct = (1:length(corr_cal)) * 343/ Fs ;
    
    idx_cal = find((ct < 1.95)) ;
    ct = ct(idx_cal) ;
    corr_cal = fft(sp_cal) ;
    corr_cal = corr_cal(idx_cal) ;
    figure; plot(ct, corr_cal);
    sp_cal = ifft(corr_cal) ;
    df_cal = Fs/length(sp_cal) ;
    f_cal = df_cal : df_cal :Fs ;
    
%% obrabotka
    [b, Fs] = wavread('2018.9.15/ogurets_test3/chan1-2.wav') ;

    [m m_idx] = max(correl(b(:,1), mls)) ;

    corr2 = correl(b(m_idx + 1: m_idx + N, 2), mls ) ;
    ct = (1:length(corr2)) * 343/ Fs ;
    figure; plot(ct, corr2) ;
    idx = (ct < 2.6) & (ct > 0.65);   
    
    corr2 = corr2(idx) ;
    corr2 = corr2(1:end-1) ;
    sp_corr2 = ifft(corr2) ./ sp_cal ;%* r(5);

    corr2 = fft(sp_corr2) ;
    ct = (1:length(corr2)) * 343/ Fs + 0.65;
    figure; plot(ct, real(corr2)) ;
    
    df = Fs/length(corr2) ;
    f_array = df:df:Fs ;
    mask = build_mask(f_array, 1000, 10000, 100, 300) ;
    sp2 = ifft(corr2) .* mask(:) ;
    figure; plot(f_array, abs(sp2)); hold all; plot(f_array, mask) ;
    
    figure; plot(ct, real(fft(sp2))/length(sp2) ) ;
