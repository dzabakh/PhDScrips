%load     control_file_1 % omni_5 %% микрофон на 1 м %%             load  omni_6 %% микрофон на 5  cм

[MLS_sample, Fs] = wavread('MLS');
[MLS_sample_3, Fs] = wavread('MLS3');


dt_sam = 1/Fs ; 

sz_sam = length(MLS_sample);
sz_sam_3 =length(MLS_sample_3);



load calibration_dir\calib_blue_red    

dt_rec = (chan1x(end)-chan1x(1)) /  (length(chan1x)-1) ;

% chanx - time, chany - singnal 
t_start = find(chan2y<-0.05,1,'first')*dt_rec ;


t_array = (t_start+(sz_sam)*dt_sam:dt_sam:t_start+(2*sz_sam-1)*dt_sam);

%%%%%%%%%%%%%%%% interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig_1 =interp1(chan1x , chan1y , t_array) ;  %%%%%%% pressure in -b
sig_2 =interp1(chan2x , chan2y , t_array) ;  %%%%%%% pressure in -a

%%%%%%%%%%%%%%%% correlation processing %%%%%%%%%%%%%%%%%%%%%%%%%

sp_1 = ifft(sig_1) ; 
sp_2 = ifft(sig_2) ;

sp_mls = ifft(MLS_sample') ;


corr_1 = real(fft(sp_1 .* conj(sp_mls))) ;
corr_2 = real(fft(sp_2 .* conj(sp_mls))) ;

%  figure;
%  plot(corr_3)
%%%%%%%%%%%%%%%% rough time mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz = sz_sam ; % 1700 ; %% size of the window (sz*10e-5*340)

corr_1 = [corr_1(end-199:end) , corr_1(1:sz-200)] ; 
corr_2 = [corr_2(end-199:end) , corr_2(1:sz-200)] ; 


%%%%%%%%%%%%%%%% filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_array = (0:sz_sam-1)*Fs /sz_sam ; 

f_high = 4000 ; 
f_width = 1500 ; 

f_max = Fs ; 

mask_1 = 0.5*(1-tanh((f_array - f_high)/f_width)) + 0.5*(1+tanh((-f_max + f_array + f_high)/f_width)) ; 

f_low = 500 ; 
f_width = 70 ; 

mask_2 = (0.5*(1+tanh((f_array - f_low)/f_width))) .*( 0.5*(1-tanh((-f_max + f_array + f_low)/f_width))) ; 

mask = mask_1 .* mask_2 ; 

corr_3_red = real(fft(ifft(corr_3) .* mask)); 

% time_array = (1:sz_sam) ; 
% 
% time_mask = exp( - time_array.^2 / 4000^2 ) ;  
% 
% corr_3_red = corr_3_red .* time_mask ;

%%%%%%%%%%%%%%%%% compressing %%%%%%%%%%%%%%

c = 340;
f_max = Fs ; 
f_array = (0:sz_sam-1)*f_max /sz_sam ; 
R_tube = 0.01875; % radius of the tube 
rho = 1.29 ; 
b = -0.05;
a = -0.03;
k_source = 2*pi*f_array/c;
sp_corr_1 = ifft(corr_1) ; 
sp_corr_2 = ifft(corr_2) ;
sp_corr_3 = ifft(corr_3_red) ; 
%%%%% two microphones

A = (-sp_corr_1.*exp(1i*k_source*b) + sp_corr_2.*exp(1i*k_source*a))./(exp(2i*k_source*a)-exp(2i*k_source*b));

W = k_source.*A*2i*pi*R_tube^2/(rho);
W(1) = W(2) ; 

sp_corr_3 = (1/sz)*(4*pi/rho)*sp_corr_3 ./ W; 

sig_3_red = 2*real(fft(sp_corr_3)) ; 

plot((1:sz_sam)*c/Fs , -sig_3_red, 'k')





% 
% 
% 
% 
% %mask = mask_1 ;
%  %mask = ones(1,3000);
% %  plot(mask)
% %%%%%%%%%%%%%%%%%%%%%%Calibration%%%%%%%%%%%%%%%%%%%
% 
%     testarray    =  zeros(1,sz);
%     testarray(1) = 1; 
%     testarray_freq = fft(testarray) ; 
%     testarray  = real(ifft(mask.*testarray_freq)) ;
%     maxval = testarray(1) ;
%     mask = mask/maxval;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% 
% sp_corr_3 = ifft(sig_3_red) ; 
% 
% sig_3_red = real(fft(sp_corr_3 .* mask)); 
% 
% 




