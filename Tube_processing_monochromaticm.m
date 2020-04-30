
clear all

ind = 0; 
load('cal_curve') ;

freqs = [500:50:3000];%[500:50:1500,1600:100:2100];
refl_coef = 0 * freqs;

spcal_interp2 = interp1(f_array, sp_2cal, freqs) ;
spcal_interp4 = interp1(f_array, sp_4cal, freqs) ;

nu = 1.56e-5 ;
r = 0.006 ;

for freq = freqs
    ind = ind + 1;
    
    folder_a = strcat('D:/ImpedanceTube/Processing/2018.4.6/mund_mono/',num2str(freq),'Hz/') ;

    chan_a = wavread([folder_a, 'chan1-2.wav']);
    chan_b = wavread([folder_a, 'chan3-4.wav']);

    sample_a = chan_a(10000+1: 110000,2); 
    etalon_a_cos = chan_a(10000+1: 110000,1) ; 
    etalon_a_sin = chan_a(10000+1-floor(48000/freq/4): 110000-floor(48000/freq/4),1) ; 

    sample_b = chan_b(10000+1: 110000,2) ;

    etalon_b_cos = chan_b(10000+1: 110000,1) ; 
    etalon_b_sin = chan_b(10000+1-floor(48000/freq/4): 110000-floor(48000/freq/4),1) ; 


    cos_a_coef = sum(sample_a .* etalon_a_cos) / sqrt(sum(etalon_a_cos.^2 + etalon_a_sin.^2)) ; 
    sin_a_coef = sum(sample_a .* etalon_a_sin) / sqrt(sum(etalon_a_cos.^2 + etalon_a_sin.^2)) ; 

    cos_b_coef = sum(sample_b .* etalon_b_cos) / sqrt(sum(etalon_b_cos.^2 + etalon_b_sin.^2)) ; 
    sin_b_coef = sum(sample_b .* etalon_b_sin) / sqrt(sum(etalon_b_cos.^2 + etalon_b_sin.^2)) ; 
    
    sp_corr_a = (cos_a_coef + 1i * sin_a_coef)/spcal_interp2(ind) ; 
    sp_corr_b = (cos_b_coef + 1i * sin_b_coef)/spcal_interp4(ind) ; 
    
    k_source = 2*pi*freq / 343 * sqrt(1 + 1i*sqrt(2)/r*exp(-1i*pi/4)*sqrt(2*nu/2/pi/freq)); 
%     a = 0.122 ; % closed end
%     b = 0.082 ; % closed end
%     a = 0.13 ;  % open end
%     b = 0.09 ; % open end

a = 0.22 ;
b = 0.18 ;

%     A = (sp_corr_a .* exp( 1i*k_source*a) - sp_corr_b .* exp( 1i*k_source*b)) ./ (exp( 2i*k_source*a) - exp( 2i*k_source*b)) ;
%     B = (sp_corr_a .* exp(-1i*k_source*a) - sp_corr_b .* exp(-1i*k_source*b)) ./ (exp(-2i*k_source*a) - exp(-2i*k_source*b)) ;     

    A = (sp_corr_b .* exp( 1i*k_source*b) - sp_corr_a .* exp( 1i*k_source*a)) ./ (exp( 2i*k_source*b) - exp( 2i*k_source*a)) ;
    B = (sp_corr_a .* exp( 1i*k_source*a + 2i*k_source*b) - sp_corr_b .* exp(1i*k_source*b + 2i*k_source*a)) ...
        ./ (exp(2i*k_source*b) - exp(2i*k_source*a)) ;

    refl_coeff(ind) = A/B ;
end
figure; hold all ;
plot(freqs, real(refl_coeff), '*-')
plot(freqs, imag(refl_coeff), '*-')
plot(freqs, abs(refl_coeff), '*-')