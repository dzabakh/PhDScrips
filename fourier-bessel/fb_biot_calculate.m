function [ Refl_coeff ] = biot_calculate( k_t, omega )

Fs = 48000 ;
h = 0.049; % layer thickness
rho_1 = 780;%1570 ;%780% granite density
rho_0 = 1.18 ;% air density
sigma = 13100;%12000 ;% flow resistivity of melamine foam (Internet)
F = 0.99 ;% porosity of melamine foam (Internet)
alpha = 1.01 ;% tortuosity of melamine foam (Internet)
nu = 0.276 ;% Poisson ratio for melamine foam (Internet)
K_a = 1.42e5 ;% bulk modulus of the fluid in the pores (air)
Lambda = 200e-6 ;%400e-6 ;   % characteristic viscous dimension
Lambda_s = 445e-6 ;%100e-6 ; % characteristic thermal dimension
eta = 18e-6 ;%air dynamic viscosity 
N = 86000*(1+0.05i) ;%shear modulus of melamine foam (Internet)
c0 = 343 ;
gamma = 1.4 ;% adiabatic coefficient for air
c_p = 1006 ;% constant pressure
khi = 26.2e-3 ;% temperature conductivity for air
k_0_s = 1.5e-9 ; %thermal permeability
K_s = 1e8;

k0 = omega / c0 ;


Pr = eta * c_p/khi; % Prandtl's number for air

M_s = 8*k_0_s/F/Lambda_s^2 ;
t_omega_s = rho_0*Pr*omega*k_0_s/eta/F ;
beta1 = gamma - (gamma-1).*(1+1./(-1i*t_omega_s).*sqrt(1-M_s/2*1i.*t_omega_s)).^(-1) ;
K_b = 2/3*N*(1+nu)/(1-2*nu) ;
K_f = K_a./beta1 ;




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


Delta1 = (P.*rho_22 + R.*rho_11 - 2*Q.*rho_12).^2 - 4*(P.*R-Q.^2).*(rho_11.*rho_22 - rho_12.^2) ;
delta_1 = sqrt(omega.^2/2./(P.*R-Q.^2).*(P.*rho_22 + R.*rho_11 - 2*Q.*rho_12 - sqrt(Delta1))) ;
delta_2 = sqrt(omega.^2/2./(P.*R-Q.^2).*(P.*rho_22 + R.*rho_11 - 2*Q.*rho_12 + sqrt(Delta1))) ;
delta_3 = sqrt( omega.^2/N .* (rho_11.*rho_22-rho_12.^2)./rho_22 ) ;

k_13 = sqrt(delta_1.^2 - k_t.^2) ;
k_23 = sqrt(delta_2.^2 - k_t.^2) ;
k_33 = sqrt(delta_3.^2 - k_t.^2) ;

mu_1 = (P.*delta_1.^2 - omega.^2.*rho_11)./(omega.^2.*rho_12 - Q.*delta_1.^2) ;
mu_2 = (P.*delta_2.^2 - omega.^2.*rho_11)./(omega.^2.*rho_12 - Q.*delta_2.^2) ;
mu_3 = -rho_11./rho_22 ;


alpha_1 = k_13 ;
alpha_2 = k_23 ;
alpha_3 = k_33 ;

beta = k_t ;

C_1 = (P + Q.*mu_1).*(beta.^2 + alpha_1.^2) - 2*N.*beta.^2 ;
C_2 = (P + Q.*mu_2).*(beta.^2 + alpha_2.^2) - 2*N.*beta.^2 ;

D_1 = (R .* mu_1 + Q) .* (beta.^2 + alpha_1.^2) ;
D_2 = (R .* mu_2 + Q) .* (beta.^2 + alpha_2.^2) ;

p_1 = cos(k_13 * h) ;
p_2 = cos(k_23 * h) ;
p_3 = cos(k_33 * h) ;

q_1 = sin(k_13 * h) ;
q_2 = sin(k_23 * h) ;
q_3 = sin(k_33 * h) ;

Delta = D_1.*(2*N.*beta.^2 + C_2) - D_2.*(2*N.*beta.^2 + C_1) ;

T_11 = (2*N.*beta.^2.*(p_2.*D_1 - p_1.*D_2) - ...
    (p_3.*(C_1.*D_2 - C_2.*D_1)))./Delta ;

T_12_top =  alpha_2.*q_1.*(mu_2.*(alpha_3.^2 - beta.^2) + ...
    2*beta.^2.*mu_3) - alpha_1.*q_2.*(mu_1.*(alpha_3.^2 - beta.^2) + ...
    2*beta.^2.*mu_3) + 2*alpha_3.*q_3.*alpha_1.*alpha_2.*(mu_1 - mu_2) ;

T_12_bottom = alpha_1 .* alpha_2 .* (mu_1 - mu_2) .* (beta.^2 + alpha_3.^2) ;

T_12 = 1i * beta .* T_12_top ./ T_12_bottom ;

T_13 = 1i * beta .* (alpha_1.*q_2 - alpha_2.*q_1)./ ...
    (alpha_1 .* alpha_2 .* (mu_1 - mu_2)) ;

T_14 = beta.*omega.*(p_1.*D_2 - p_2.*D_1 - p_3.*(D_2 - D_1))./Delta ;

T_15 = 1i*omega./N./(beta.^2 + alpha_3.^2) .* ...
    (beta.^2.*q_1.*(mu_2-mu_3)./alpha_1./(mu_2 - mu_1) +...
    beta.^2.*q_2.*(mu_1 - mu_3)./alpha_2./(mu_1 - mu_2) + alpha_3.*q_3) ;

T_16 = beta .* omega.* (p_2.*C_1 - p_1.*C_2 - p_3.*(C_1 - C_2) +...
    2*N.*beta.^2.*(p_2 - p_1))./Delta ;

T_21 = 1i*beta./Delta.*(2*N.*(alpha_1.*q_1.*D_2 - alpha_2.*q_2.*D_1) - ...
    q_3./alpha_3.*(C_1.*D_2 - C_2.*D_1)) ;

T_22 = (p_2.*(mu_1.*(alpha_3.^2 - beta.^2) + 2*beta.^2.*mu_3) - ...
        p_1.*(mu_2.*(alpha_3.^2 - beta.^2) + 2*beta.^2.*mu_3) + ...
        2*beta.^2.*p_3.*(mu_1 - mu_2))./...
        ((mu_1 - mu_2).*(beta.^2 + alpha_3.^2)) ;
    
T_23 = (p_1 - p_3)./(mu_1 - mu_2) ;

T_24 = -1i.*omega./Delta .* (alpha_1.*q_1.*D_2 - alpha_2.*q_2.*D_1 + ...
    beta.^2.*q_3./alpha_3.*(D_2 - D_1)) ;

T_25 = - omega.*beta./N./(beta.^2 + alpha_3.^2).*(p_1.*(mu_2 - mu_3)./...
    (mu_2 - mu_1) + p_2.*(mu_1 - mu_3)./(mu_1 - mu_2) - p_3) ;

T_26 = 1i*omega./Delta.*(alpha_1.*q_1.*(C_2 + 2*N.*beta.^2) - ...
    alpha_2.*q_2.*(C_1 + 2*N.*beta.^2)-...
    q_3.*beta.^2/alpha_3.*(C_1 - C_2)) ;

T_31 = 1i*beta./Delta.*(2*N.*(alpha_1.*mu_1.*q_1.*D_2 - ...
    alpha_2.*mu_2.*q_2.*D_1) - mu_3.*q_3./alpha_3.*(C_1.*D_2 - C_2.*D_1)) ;

T_32 = (-mu_1.*p_1.*(mu_2.*(alpha_3.^2 - beta.^2) + 2.*beta.^2.*mu_3) + ...
    mu_2.*p_2.*(mu_1.*(alpha_3.^2 - beta.^2) + 2*beta.^2.*mu_3) + ...
    2*beta.^2.*mu_3.*p_3.*(mu_1 - mu_2)) ./...
    ((mu_1 - mu_2).*(alpha_3.^2 + beta.^2)) ;

T_33 = (mu_1.*p_1 - mu_2.*p_2)./(mu_1 - mu_2) ;

T_34 = 1i*omega./Delta.*(-alpha_1.*mu_1.*q_1.*D_2 + alpha_2.*mu_2.*q_2.*D_1 + ...
    beta.^2.*mu_3.*q_3./alpha_3.*(D_1 - D_2)) ;

T_35 = -beta.*omega/N./(alpha_3.^2 + beta.^2).*(p_1.*mu_1.*(mu_2 - mu_3)./(mu_2 - mu_1) + ...
    p_2.*mu_2.*(mu_1 - mu_3)./(mu_1 - mu_2) - p_3.*mu_3) ;

T_36 = 1i*omega./Delta.*(mu_1.*alpha_1.*q_1.*(C_2 + 2*N*beta.^2) - ...
    mu_2.*alpha_2.*q_2.*(C_1 + 2*N*beta.^2) - ...
    beta.^2/alpha_3.*mu_3.*q_3.*(C_1 - C_2)) ;

T_41 = 2*N*beta./omega./Delta.*(C_1.*p_1.*D_2 - C_2.*p_2.*D_1 - ...
    p_3.*(C_1.*D_2 - C_2.*D_1)) ;

T_42_top = C_1.*q_1.*alpha_2.*(mu_2.*(alpha_3.^2 - beta.^2) + 2*beta.^2.*mu_3) - ...
    C_2.*q_2.*alpha_1.*(mu_1.*(alpha_3.^2 - beta.^2) + 2*beta.^2.*mu_3) - ...
    4*N*alpha_3.*beta.^2.*alpha_1.*alpha_2.*(mu_1 - mu_2).*q_3 ;

T_42_bottom = alpha_1.*alpha_2.*omega.*(beta.^2 + alpha_3.^2).*(mu_1 - mu_2) ;

T_42 = -1i*T_42_top./T_42_bottom ;

T_43 = 1i*(alpha_2.*C_1.*q_1 - alpha_1.*C_2.*q_2)./...
    (omega.*alpha_1.*alpha_2.*(mu_1 - mu_2)) ;

T_44 = (-p_1.*C_1.*D_2 + p_2.*C_2.*D_1 - 2*N.*beta.^2.*p_3.*(D_2 - D_1))...
    ./Delta ;

T_45 = -1i*beta./(beta.^2 + alpha_3.^2).*...
    (C_1.*q_1.*(mu_2 - mu_3)./N./alpha_1./(mu_2 - mu_1) + ...
    C_2.*q_2.*(mu_1 - mu_3)./N./alpha_2./(mu_1 - mu_2) - 2.*q_3.*alpha_3) ;

T_46 = (p_1.*C_1.*(C_2 + 2.*N.*beta.^2) - p_2.*C_2.*(C_1+2*N.*beta.^2) - ...
    2*N.*beta.^2.*p_3.*(C_1 - C_2)) ./ Delta ;

T_51 = 1i*N.*beta.^2/Delta./omega .* (4*N*alpha_1.*q_1.*D_2 - 4*N*alpha_2.*q_2.*D_1 ...
    - q_3.*(alpha_3.^2 - beta.^2)./beta.^2./alpha_3.*(C_1.*D_2 - C_2.*D_1)) ;

T_52_top = 2*N*beta.*p_1.* (mu_2.*(alpha_3.^2 - beta.^2) + 2*beta.^2.*mu_3) - ...
    2*N*beta.*p_2.*(mu_1.*(alpha_3.^2 - beta.^2) + 2*beta.^2.*mu_3) + ...
    2*N*beta.*p_3.*(alpha_3.^2 - beta.^2).*(mu_1 - mu_2) ;

T_52_bottom = omega.*(mu_1 - mu_2) .* (beta.^2 + alpha_3.^2) ;

T_52 = T_52_top ./ T_52_bottom ;

T_53 = -2*N*beta./omega./(mu_1 - mu_2).*(p_1 - p_2) ;

T_54 = 2*1i*N.*beta./Delta.*(alpha_1.*q_1.*D_2 - alpha_2.*q_2.*D_1 - ...
    q_3/2.*(alpha_3.^2 - beta.^2)./alpha_3.*(D_2 - D_1)) ;

T_55 = 2*beta.^2./(beta.^2 + alpha_3.^2).*(p_1.*(mu_2 - mu_3)./(mu_2 - mu_1) + ...
    p_2.*(mu_1 - mu_3)./(mu_1 - mu_2) + p_3.*(alpha_3.^2 - beta.^2)/2./(beta.^2)) ;

T_56 = -2i*N*beta./Delta.*(alpha_1.*q_1.*(C_2 + 2*N.*beta.^2) - alpha_2.*q_2.*(C_1 ...
    + 2*N*beta.^2) + q_3/2.*(alpha_3.^2 - beta.^2)./alpha_3.*(C_1 - C_2)) ;

T_61 = 2*N*beta.*D_1.*D_2./omega./Delta.*(p_1 - p_2) ;

T_62_top = alpha_2.*q_1.*D_1.*(mu_2.*(alpha_3.^2 - beta.^2) + 2*beta.^2.*mu_3) - ...
    alpha_1.*q_2.*D_2.*(mu_1.*(alpha_3.^2 - beta.^2) + 2*beta.^2.*mu_3);

T_62_bottom = alpha_1.*alpha_2 .* (mu_1 - mu_2).*(beta.^2 + alpha_3.^2) ;

T_62 = -1i./omega .* T_62_top ./ T_62_bottom ;

T_63 = 1i ./ omega./(mu_1 - mu_2).*(q_1.*D_1./alpha_1 - q_2.*D_2./alpha_2) ;

T_64 = - D_1.*D_2./Delta.* (p_1 - p_2) ;

T_65 = -1i*beta./N./(beta.^2 + alpha_3.^2).*(q_1.*D_1./alpha_1.*(mu_2 - mu_3)./...
    (mu_2 - mu_1) + q_2.*D_2./alpha_2.*(mu_1 - mu_3)./(mu_1 - mu_2)) ;

T_66 = (p_1.*D_1.*(C_2 + 2*N*beta.^2) - p_2.*D_2.*(C_1 + 2*N*beta.^2)) ./Delta ;


Delta_1 = T_24.*T_55 - T_25.*T_54 ;
Delta_2 = T_26.*T_55 - T_25.*T_56 ;
Delta_3 = T_34.*T_55 - T_35.*T_54 ;
Delta_4 = T_36.*T_55 - T_35.*T_56 ;
Delta_5 = T_44.*T_55 - T_45.*T_54 ;
Delta_6 = T_46.*T_55 - T_45.*T_56 ;
Delta_7 = T_64.*T_55 - T_65.*T_54 ;
Delta_8 = T_66.*T_55 - T_65.*T_56 ;

D = ((1-F).*Delta_1 + F.*Delta_3).*((1-F).*Delta_8 - F.*Delta_6) - ...
    ((1-F).*Delta_2 + F.*Delta_4).*((1-F).*Delta_7 - F.*Delta_5) ;


Z = (Delta_6.*Delta_7 - Delta_5.*Delta_8) ./ D ;

gamma_0 = sqrt(k0.^2 - k_t.^2) ;
Refl_coeff = (1 - rho_0*omega./Z./gamma_0)./(1 + rho_0*omega./Z./gamma_0) ;


end

