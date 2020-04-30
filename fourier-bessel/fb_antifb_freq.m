clear all;
load('M_matrix2');
load('M_matrix');

load('melamin\R0_mel') ;
load('melamin\R20_mel') ;
load('melamin\R30_mel') ;
load('melamin\R40_mel') ;
load('melamin\R45_mel') ;
load('melamin\R50_mel') ;
load('melamin\R60_mel') ;
load('melamin\R70_mel') ;
R90_mel = R60_mel * 0 ;

rho0 = 1.21 ;
R_j_raw = zeros(length(f_array), 5);%length(theta_exp)-1) ; % \tilde{R(\phi_i)}
R_j = zeros(length(f_array), 5);%length(theta_exp)-1) ;     % \tilde{R(\phi_i)}


for f = 50:626
      omega = 2*pi*f_array(f) ;
      u_n_mel  = [R0_mel(f);  R20_mel(f); R30_mel(f); R45_mel(f); R60_mel(f)] ;
%       u_n_wall = [R0_wall(f); R10_wall(f); R20_wall(f); R30_wall(f); R45_wall(f); R60_wall(f)] ;
      R_j_raw(f,:) = pinv(M(:,:,f)) * u_n_mel * (rho0/(4*pi))* sz;
  end 

%% Masking and norming
f_high = 4500 ; 
f_width =  500 ;
f_s = f_array(end);
    
mask = 0.5*(1-tanh((f_array - f_high)/f_width)) + 0.5*(1+tanh((-f_s + f_array + f_high)/f_width)) ; 

mask(length(mask)/2 + 1: end) = 0; 

    
%%

for i = 1:5
   R_j(:,i) = R_j_raw(:,i) .* mask.' ;
end    

figure;
plot(f_array, abs(R_j(:,1)));

% figure;
% plot(f_array, R_j(:,2));
% 
% figure;
% plot(f_array, R_j(:,3));
% 
% figure;
% plot(f_array, R_j(:,4));
% 
% figure;
% plot(f_array, R_j(:,5));
% 
% figure;
% plot(f_array, R_j(:,6));
% 
% figure;
% plot(f_array, R_j(:,7));
% 
% figure;
% plot(f_array, R_j(:,8));

