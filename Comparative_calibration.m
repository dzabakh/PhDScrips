clear all
[MLS_sample, Fs] = wavread('MLS');
load('cal_coeffs') ;

    sz_sam = length(MLS_sample);

%     [a21, Fs_rec] = wavread('2015.12.14/comp_calibr_bk_a2/chan5-6.wav') ;
%     [bk, Fs_rec] = wavread('2015.12.14/comp_calibr_bk_a2/chan3-4.wav') ;
%         
%     [a22, Fs_rec] = wavread('2015.12.14/comp_calibr_a1_a2/chan5-6.wav') ;
%     [a1, Fs_rec] = wavread('2015.12.14/comp_calibr_a1_a2/chan7-8.wav') ;
    
    [a21, Fs_rec] = wavread('2015.12.22/comp_calibr2/chan5-6.wav') ;
    [bk, Fs_rec] = wavread('2015.12.22/comp_calibr2/chan3-4.wav') ;
        
    [a22, Fs_rec] = wavread('2015.12.22/comp_calibr2/chan5-6.wav') ;
    [a1, Fs_rec] = wavread('2015.12.22/comp_calibr2/chan7-8.wav') ;
        


    sz_sam = length(MLS_sample);
    
    xcorr_a21 = fft(ifft(a21(1:sz_sam,2)) .* conj(ifft(MLS_sample))) ;
    [maxval_a21, t_start_a21] = max(abs(xcorr_a21)) ;
    t_start_a21 = t_start_a21 - 10 ;
    t_end_a21 = t_start_a21 + sz_sam - 1 ;
    
    xcorr_bk = fft(ifft(bk(1:sz_sam,2)) .* conj(ifft(MLS_sample))) ;
    [maxval_bk, t_start_bk] = max(abs(xcorr_bk)) ;
    t_start_bk = t_start_bk - 10 ;
    t_end_bk = t_start_bk + sz_sam - 1 ;
    
    xcorr_a22 = fft(ifft(a22(1:sz_sam,2)) .* conj(ifft(MLS_sample))) ;
    [maxval_a22, t_start_a22] = max(abs(xcorr_a22)) ;
    t_start_a22 = t_start_a22 - 10 ;
    t_end_a22 = t_start_a22 + sz_sam - 1 ;
    
    xcorr_a1 = fft(ifft(a1(1:sz_sam,2)) .* conj(ifft(MLS_sample))) ;
    [maxval_a1, t_start_a1] = max(abs(xcorr_a1)) ;
    t_start_a1 = t_start_a1 - 10 ;
    t_end_a1 = t_start_a1 + sz_sam - 1 ;
    
    bk = bk(t_start_bk : t_end_bk, 2)' ;
    a1 = a1(t_start_a1 : t_end_a1, 2)' ;
    a21 = a21(t_start_a21 : t_end_a21, 2)' ;
    a22 = a22(t_start_a22 : t_end_a22, 2)' ;
    
    bk  = bk  * coeff_bk ;
    a1  = a1  * coeff1   ;
    a21 = a21 * coeff2   ;
    a22 = a22 * coeff2   ;
    
    sp_bk = ifft(bk) ; 
    sp_21 = ifft(a21);  
    sp_1  = ifft(a1) ;   
    sp_22 = ifft(a22);


    sp_mls = ifft(MLS_sample') ; 

    corr_bk = real(fft(sp_bk .* conj(sp_mls))) ;
    corr_21 = real(fft(sp_21 .* conj(sp_mls))) ;
    corr_1  = real(fft(sp_1  .* conj(sp_mls))) ;
    corr_22 = real(fft(sp_22 .* conj(sp_mls))) ;
  
    sz = 10000 ;
    corr_bk = corr_bk(1:sz) ; 
    corr_21 = corr_21(1:sz) ;
    corr_1  = corr_1(1:sz)  ;
    corr_22 = corr_22(1:sz) ;
    
    mask = 0.5*(1-tanh(((1:sz)-500)/100));
    
    corr_bk = corr_bk .* mask ; 
    corr_21 = corr_21 .* mask ; 
    corr_1  = corr_1  .* mask ; 
    corr_22 = corr_22 .* mask ; 
    
    c = 340;
    f_max = Fs ; 
    f_array = (0:sz-1)*f_max /sz ; 
    
    sp_bk = ifft(corr_bk) ; 
    sp_21 = ifft(corr_21)  ; 
    sp_1  = ifft(corr_1)   ; 
    sp_22 = ifft(corr_22)  ; 
    
    figure
    hold all
    plot(f_array, real(sp_21./sp_bk)) ;
    plot(f_array, imag(sp_21./sp_bk)) ;
    
    figure
    hold all
    plot(f_array, real(sp_22./sp_1)) ;
    plot(f_array, imag(sp_22./sp_1)) ;
    
    cal_f_array = f_array ; 
    cal_curve_bk_a2 = sp_21 ./ sp_bk ; 
    cal_curve_a1_a2 = sp_22  ./ sp_1 ; 
        
    figure ;
    hold all;
    plot(real(fft(sp_21) )) ;
    plot(real(fft(sp_bk .* cal_curve_bk_a2))) ;
    
    save('cal_curves.mat', 'cal_f_array', 'cal_curve_a1_a2', 'cal_curve_bk_a2') ;
    