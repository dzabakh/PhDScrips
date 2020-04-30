clear all;

sz = 3000 ;
c0 = 343 ;

freq = (1:sz-1)*48000 /sz ; 

omega = 2*pi*freq ;
k0 = omega / c0 ;
theta = 1/180*pi;%(1.001:1:90.001) / 180 * pi; 
kt1 = k0 * sin(theta) ;
% kt2 = [(k0+k0/1000):k0/1000:(k0+k0/100),(k0+2*k0/100):k0/100:(4*k0)]; 
% k_t = [kt1,kt2];
k_t = kt1;

% k_t = 0.48;
refl_coeff = k_t*0 ;
for i = 1:length(k_t)
    refl_coeff(i) = biot_calculate(k_t(i), omega(i)) ;
end
% figure; plot(freq, real(refl_coeff)) ;
figure; plot(freq, real(refl_coeff)) ;