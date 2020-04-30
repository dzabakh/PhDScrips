% clear all;

Fs = 48000 ;
h = 0.05; % layer thickness
rho_1 = 780 ;% granite density
rho_0 = 1.21 ;% air density
sigma = 12000 ;% flow resistivity of melamine foam (Internet)
F = 0.99 ;% porosity of melamine foam (Internet)
alpha = 1.01 ;% tortuosity of melamine foam (Internet)
nu = 0.276 ;% Poisson ratio for melamine foam (Internet)
K_a = 1.42e5 ;% bulk modulus of the fluid in the pores (air)
Lambda = 100e-6 ;% characteristic viscous dimension
Lambda_s = 400e-6 ; % characteristic thermal dimension
eta = 15e-6 ;%air dynamic viscosity 
N = 86000*(1+0.05i) ;%shear modulus of melamine foam (Internet)
c0 = 343 ;
gamma = 1.4 ;% adiabatic coefficient for air
c_p = 1006 ;% constant pressure
khi = 26.2e-3 ;% temperature conductivity for air
k_0_s = 1.5e-9 ; %thermal permeability
K_s = 1e8;

sz = 5000 ;


freq = 2000 /sz ; 


omega = 2*pi*freq ;
theta = asin(0:0.001:10) ;%(1:1:90) / 180 * pi; %0/180*pi;


Pr = eta * c_p/khi; % Prandtl's number for air

M_s = 8*k_0_s/F/Lambda_s^2 ;
t_omega_s = rho_0*Pr*omega*k_0_s/eta/F ;
beta = gamma - (gamma-1)*(1+1./(-1i*t_omega_s).*sqrt(1-M_s/2*1i.*t_omega_s)).^(-1) ;
K_b = 2/3*N*(1+nu)/(1-2*nu) ;
K_f = K_a./beta ;




% P, Q, R for material with incompressible skeleton
P = 4/3*N + K_b + (1-F)^2/F*K_f ; %( (1-F)*(1-F-K_b/K_s)*K_s + F*K_s./K_f*K_b)./(1-F-K_b/K_s+F*K_s./K_f) + 4/3*N ;
Q = K_f*(1-F) ;%(1-F-K_b/K_s)*F*K_s./(1-F-K_b/K_s + F*K_s./K_f) ;
R = F*K_f ;%F^2*K_s./(1-F-K_b/K_s + F*K_s./K_f) ;

G = sqrt(1-4i*alpha^2*eta*rho_0.*omega/sigma^2/Lambda^2/F^2);
b = -sigma./omega*F^2.*G ;

rho_a = F*rho_0*(alpha - 1);
rho_11 = rho_1 + rho_a - 1i*b ;
rho_12 = -rho_a + 1i*b ;
rho_22 = F*rho_0 + rho_a - 1i*b ;





Delta = (P.*rho_22 + R.*rho_11 - 2*Q.*rho_12).^2 - 4*(P.*R-Q.^2).*(rho_11.*rho_22 - rho_12.^2) ;
delta_1 = sqrt(omega.^2/2./(P.*R-Q.^2).*(P.*rho_22 + R.*rho_11 - 2*Q.*rho_12 - sqrt(Delta))) ;
delta_2 = sqrt(omega.^2/2./(P.*R-Q.^2).*(P.*rho_22 + R.*rho_11 - 2*Q.*rho_12 + sqrt(Delta))) ;
delta_3 = sqrt( omega.^2/N .* (rho_11.*rho_22-rho_12.^2)./rho_22 ) ;

mu_1 = (P.*delta_1.^2 - omega.^2.*rho_11)./(omega.^2.*rho_12 - Q.*delta_1.^2) ;
mu_2 = (P.*delta_2.^2 - omega.^2.*rho_11)./(omega.^2.*rho_12 - Q.*delta_2.^2) ;
mu_3 = -rho_11./rho_22 ;

Z_1f = (R + Q./mu_1).*delta_1/F./omega ;
Z_2f = (R + Q./mu_2).*delta_2/F./omega ;

Z_1s = (P + Q.*mu_1).*delta_1./omega ;
Z_2s = (P + Q.*mu_2).*delta_2./omega ;

D = (1-F+F*mu_2) .* (Z_1s - (1-F)*Z_1f.*mu_1).*tan(delta_2*h) + (1-F+F*mu_1).*(Z_2f.*mu_2*(1-F) - Z_2s).*tan(delta_1*h) ;
Z = 1i*(Z_1s.*Z_2f.*mu_2 - Z_2s.*Z_1f.*mu_1)./D ;%%%%%%%%%??????????

%  figure;
% 
%  plot(freq,real(Z)/rho_0/c0) ; 
%  hold all;
%  plot(freq,imag(Z)/rho_0/c0) ;





k0 = omega/c0 ;
kii = k0 .* sin(theta) ;
theta = asin(kii/k0) ;
gamma_0 = sqrt(k0.^2 - kii.^2) ;
R = (1 - rho_0*omega./Z./gamma_0)./(1 + rho_0*omega./Z./gamma_0) ;


R = R(2:end);


%%%%%%%%%%FILTER%%%%%%%%%%%%%%%%%

 
c = 340 ;

f_max = Fs ; 
f_array = freq ;%(0:sz-1)*f_max /sz ; 
df = f_max/sz ;
k_source = 2*pi*f_array/c;

f_high = 4500 ; 
f_width =  500 ; 
mask_1 = 0.5*(1-tanh((f_array - f_high)/f_width)) + 0.5*(1+tanh((-f_max + f_array + f_high)/f_width)) ; 

    
f_low = 500 ; 
f_width = 100 ; 
mask_2 = (0.5*(1+tanh((f_array - f_low)/f_width))) .*( 0.5*(1-tanh((-f_max + f_array + f_low)/f_width))) ; 

mask = mask_1;%.*mask_2 ; 


 R = abs(R);%.*mask(2:end) ;

figure;
plot(kii(2:end)/k0, R, 'LineWidth',2,'Color','red') ;
% plot(theta(1:end-1)/pi*180, R, 'LineWidth',2,'Color','red') ;
% theta(1:end-1)/pi*180
