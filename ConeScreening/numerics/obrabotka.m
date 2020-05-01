clear all

[mls, Fs] = audioread('mls_long.wav') ;
mls = mls(:,1) ;
N = length(mls) ;    
sz = 3000;

%% kalibrovka
    cali_file = '/Volumes/hdd/dzabakh/Dropbox/SaintGobain/ConeScreening/2018.5.15/cal50/chan7-8.wav' ;
    cal = audioread(cali_file) ;
    [m, m_idx] = max(correl(cal(:,1), mls)) ;
    corr_cal = real(correl(cal(m_idx + 1: m_idx + N, 2), mls)) ;
    
   ct = (1:sz) * 343/ Fs ;
    corr_cal(ct<0.25)=0;
    corr_cal(ct>1.2)=0;
    corr_cal = corr_cal(1:sz);
    
    df = Fs/length(corr_cal) ;
    f_array = df:df:Fs ;
    k = 2*pi*f_array/343 ;
    
    sp_cal = ifft(corr_cal) .* exp(-1i*k.'*0.5) * 0.5;
    
    
%% obrabotka
   % y = 0.1 : 0.1 : 1 ; %% shadow  
    y = -0.1 : -0.1 : -0.5 ; %%light
    x0 = 1;
    l_s = 0.8;
    alpha = 0.0456;
    r = sqrt((l_s+x0)^2+ (l_s*alpha+y).^2);
    n_exp =length(r);

    sig_int_sp = zeros(n_exp,length(corr_cal));
    for n_cur = 1:n_exp
%         cur_filename =  ['/Volumes/hdd/dzabakh/Dropbox/SaintGobain/ConeScreening/2018.5.15/tenj',num2str(n_cur),'/chan7-8.wav'];
        cur_filename =  ['/Volumes/hdd/dzabakh/Dropbox/SaintGobain/ConeScreening/2018.5.15/svet',num2str(n_cur),'/chan7-8.wav'];
        [b, Fs] = audioread(cur_filename) ;

        [m, m_idx] = max(correl(b(:,1), mls)) ;

        corr2 = real(correl(b(m_idx + 1: m_idx + N, 2), mls )) ;
        ct_temp = (1:length(corr2)) * 343/ Fs ;

        corr2(ct_temp>2.78)=0;
        corr2(ct_temp<1.6)=0;
        corr2 = corr2(1:sz);

        sp_corr2 = ifft(corr2) ./ sp_cal * r(n_cur) .* exp(-1i * k.' * r(n_cur)) ;
        sp2 = sp_corr2 ;
        sig_int_sp(n_cur,:) = sp2;
    end
    
    save('svet','sig_int_sp','r','f_array')
