clear all;

%folder = '2019.6.13\empty_room_small_speaker\';
%folder = '2019.6.19\on_box_at_room_part2\';
 folder = '2019.6.19\calibration_1m\';
%folder = '2019.6.13\twitter_calibration_10cm\';
 load('source_cal_coef_19_06')

sz = 3000 ;
c = 343;
R_mic = 0.12; % distance
rho = 1.21 ;

chan12 = audioread([folder, 'chan1-2.wav']);


[mls, Fs] = audioread('mls_long.wav');
mls = mls(:,1) ;
sz_mls = length(mls) ;

xcorr12 = real(fft(ifft(chan12(1:sz_mls,1)) .* conj(ifft(mls)))) ;


[max12, t_start12] = max(xcorr12) ;


t_end12 = t_start12 + sz_mls - 1 ;



probemic = chan12(t_start12:t_end12, 2).';


sp_probe = ifft(probemic) ;
sp_mls = ifft(mls.') ;
corr_probe = real(fft(sp_probe .* conj(sp_mls))) ;
%  figure;
% % plot(corr_a1)
% % hold all
% % plot(corr_a2)
% % hold all
% 
%  plot(0:c/Fs:c*(length(corr_probe)-1)/Fs,corr_probe)

%plot(0:c/Fs:c*(length(chan12)-1)/Fs,chan12(:, 2))
%%



f_max = Fs ; 
f_array = (0:sz-1)*f_max /sz ;
t_ar = 0:1/Fs:(sz-1)/Fs; 
df = f_max/sz ;
k_ar = f_array*2*pi/c;



corr_probe = corr_probe(1:sz) ;


    
%% time masking %%%%

% t_down = 0.5/c;
% dt_down = 0.1/c;
% mask_time =  1-1/2 * ( tanh((t_ar- t_down)/dt_down) + 1) ;
% 
%  corr_probe = mask_time.*corr_probe;


%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp_corr_probe = ifft(corr_probe);
%%%%%%%%%%%%%% masking
f_high = 18000 ; 
f_width =  2000 ; 
    
mask_1 = 0.5*(1-tanh((f_array - f_high)/f_width)) + 0.5*(1+tanh((-f_max + f_array + f_high)/f_width)) ; 

    
f_low = 3000 ; 
f_width = 500 ; 

mask_2 = (0.5*(1+tanh((f_array - f_low)/f_width))) .*( 0.5*(1-tanh((-f_max + f_array + f_low)/f_width))) ; 

mask = mask_1.* mask_2 ; 
mask(length(mask)/2 + 1: end) = 0; 


  %%%%%%%%%%%%%%%%%%%%%Calibration%%%%%%%%%%%%%%%%%%%

    %delta = 0 * mask ; 
    delta_freq = exp(1i*k_ar); 
    delta  = real(ifft(mask.*delta_freq)) ;
    normirovka = max(delta);
    mask = mask/normirovka;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    



% source_cal_coef = 1./sp_corr_probe.*exp(1i*R_mic*k_ar)/R_mic/sz;
 
 sp_corr_probe = sp_corr_probe.* mask;
 
%%   
%    figure;
%   plot(c*t_ar,corr_probe)
%   hold all
%   plot(c*t_ar, mask_time)
  figure;
  plot(c*t_ar, real(fft(sp_corr_probe.*source_cal_coef)))
 %title('$|P_3|$','Interpreter','latex','FontSize',14)
  xlabel('$ct$, m','Interpreter','latex','FontSize',14)
  ylabel('Амплитуда ','Interpreter','latex','FontSize',14)
% 
%       save('source_cal_coef_19_06','source_cal_coef','f_array')


