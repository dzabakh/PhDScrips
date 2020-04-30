% clear all;
load('M_matrix2');
load('M_matrix');

% load('R0_wall') ;
% load('R10_wall') ;
% load('R20_wall') ;
% load('R30_wall') ;
% load('R40_wall') ;
% load('R45_wall') ;
% load('R50_wall') ;
% load('R60_wall') ;
% R90_wall = R60_wall * 0 ;

load('melamin\R0_mel') ;
load('melamin\R10_mel') ;
load('melamin\R20_mel') ;
load('melamin\R30_mel') ;
load('melamin\R40_mel') ;
load('melamin\R45_mel') ;
load('melamin\R50_mel') ;
load('melamin\R60_mel') ;
load('melamin\R70_mel') ;
R90_mel = R60_mel * 0 ;

f_plot = 4000 ;

df = f_array(2) - f_array(1) ;

rho0 = 1.18 ;
% theta_exp = [0, 20, 30, 45, 60]/180*pi ;
% theta = (0:1:60) /180*pi ;
theta_exp = [1, 10, 20, 30, 40, 45, 50, 60, 70]/180*pi ;
R_j_raw = zeros(length(f_array), length(theta_exp));


for f = 105:626
%       u_n_mel  = [R0_mel(f); R20_mel(f); R30_mel(f); R45_mel(f); R60_mel(f)] ;
    u_n_mel  = [R0_mel(f); R10_mel(f); R20_mel(f); R30_mel(f); R40_mel(f); R45_mel(f); R50_mel(f); R60_mel(f); R70_mel(f)] ;
    Mm1 = pinv(M(:,:,f));

    k_0 = f_array(f)*2*pi/c0 ;
    dk = k_0/100 ;
    k_ii = ((k_0 + 2*dk) : dk: 10*k_0) ;
    theta = asin(k_ii/k_0) ;
    R_biot = biot_func(f_array(f), theta) ;
    yadro = zeros(length(theta_exp), length(k_ii)) ;
    for ind = 1:length(theta_exp)
        yadro(ind, :) = besselj(0, 2*k_ii* sin(theta_exp(ind))) .* k_ii./sqrt(k_0^2 - k_ii.^2) .* exp(1i*sqrt(k0^2 - k_ii.^2)*2*cos(theta_exp(ind))) ;
        integr2(ind) = sum(R_biot .* yadro(ind,:) * dk)   ;  
    end
    R_j_raw(f,:) = Mm1 * (u_n_mel ) * (rho0/(4*pi))* sz;
    
end 

theta_exp = [theta_exp, pi/2] ;
R_j_raw = [R_j_raw, -ones(length(f_array),1)] ;
  

f_high = 4500 ; 
f_width =  500 ;
f_s = f_array(end);
    
mask = 0.5*(1-tanh((f_array - f_high)/f_width)) + 0.5*(1+tanh((-f_s + f_array + f_high)/f_width)) ; 

mask(length(mask)/2 + 1: end) = 0; 





for j = 1:length(theta_exp)
   R_j_raw(:,j) = R_j_raw(:,j) .* mask.'; 
end    
    

% [X,Y] = meshgrid(theta_exp/pi*180, f_array(105:626)) ;
% figure; 
% surf(X, Y, abs(R_j_raw(105:626,:)));

to_plot = real(R_j_raw(floor(f_plot/df),:));
figure;
plot(theta_exp/pi*180, to_plot , '-*');
title(['f = ',num2str(f_plot)]);
ylabel('R');
