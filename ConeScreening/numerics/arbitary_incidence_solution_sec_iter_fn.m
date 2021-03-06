clear all

 
%% setting parameters
f_ar = 3e3:1e3:1e4;
u_theor_svet = 1:length(f_ar);
u_theor_tenj = 1:length(f_ar);
%f = 5000;
c = 343;
%k = 2*pi*f/c;
alpha = 0.0456;
%ka2 = k*alpha^2;
x_end = 0.8; 
N_x = 5;
n_im = 100;
h_im = 0.005;
n_pres = 1000;
N_iter = 20;
dx = x_end/N_x;
x_ar = dx:dx:x_end;
dphi = 0.05;
phi_ar = 0:dphi:(2*pi-dphi);
th_in = atan(0.1);
x_in  = -1;
rho_in = 0.1;
phi_in = 0;


%% calculating kernel
for cur_f_n = 1:length(f_ar) 
Kernel_full = zeros(length(phi_ar)-1,length(phi_ar),N_x,N_x);
k = 2*pi*f_ar(cur_f_n)/c;

wb = waitbar(0 , 'smoking bamboo') ;
for cur_n_1=1:length(x_ar)
    x = x_ar(cur_n_1);
    waitbar((cur_n_1-1)/(length(x_ar)-1))
    for cur_n_phi = 1:length(phi_ar)
        phi = phi_ar(cur_n_phi);
        for cur_n_2 = 1:cur_n_1                      
            if (cur_n_2==1)
                x3 = x_ar(cur_n_2+1);
                x2 = x_ar(cur_n_2);
                Kernel_full(:,cur_n_phi,cur_n_1,cur_n_2) = Kij_full(k,x,x3,x2,phi_ar,phi,alpha,n_pres,n_im,h_im,1);
            elseif (cur_n_2 == cur_n_1)%&&(cur_n_2~=1)
                x2 = x_ar(cur_n_2);
                x1 = x_ar(cur_n_2-1);
                Kernel_full(:,cur_n_phi,cur_n_1,cur_n_2) = Kij_full(k,x,x2,x1,phi_ar,phi,alpha,n_pres,n_im,h_im,0);
            else
                x3 = x_ar(cur_n_2+1);
                x2 = x_ar(cur_n_2);
                x1 = x_ar(cur_n_2-1);
                Kernel_full(:,cur_n_phi,cur_n_1,cur_n_2) = Kij_full(k,x,x2,x1,phi_ar,phi,alpha,n_pres,n_im,h_im,0)...
                + Kij_full(k,x,x3,x2,phi_ar,phi,alpha,n_pres,n_im,h_im,1);       
            end
        end
    end
end
close(wb)
Ker_add = Kernel_full(:,:,2:end,1);
Kernel_full = Kernel_full(:,:,2:end,2:end);
% save('Kernel_full_1','Ker_add','Kernel_full')
%%  making iterations
% load('Kernel_full_1')
% U_in_full = ones(length(phi_ar)-1,length(phi_ar)-1,N_x-1,N_x-1);
% U_in_1 = ones(length(phi_ar)-1,N_x-1);
% U_in_2 = ones(length(phi_ar),N_x-1);


 
ones_phi = ones(length(phi_ar),1);
x_grid = ones_phi*x_ar(2:end);
rho_grid = (x_grid*alpha);
phi_in_grid = (phi_ar.'-phi_in)*ones(1,length(x_ar)-1);
delta_rho_grid2 = rho_grid.^2 + rho_in^2 - 2*rho_grid*rho_in.*cos(phi_in_grid);
    
% U_in_2= exp(-1i*k*x_grid.*(th_in^2/2+th_in*alpha.*cos(phi_in_grid))); % plane wave

U_in_2 = k./(2i*pi*(x_grid-x_in)).*exp(0.5i*k.*delta_rho_grid2./(x_grid-x_in)); % point source

U_in_1 = (U_in_2(2:end,:)+U_in_2(1:end-1,:))/2;
u = 0*U_in_2;
u_prev = U_in_1;
%  figure; 
       for n_cur_it = 1:N_iter
           for cur_n_x = 1:(length(x_ar)-1)
               for cur_n_phi = 1:(length(phi_ar)) 
                 u(cur_n_phi,cur_n_x) = sum(sum((squeeze(Kernel_full(:,cur_n_phi,cur_n_x,:)).*u_prev),1),2) + squeeze(sum(Ker_add(:,cur_n_phi,cur_n_x,:),1))+2.*U_in_2(cur_n_phi,cur_n_x);
               end
           end
           u_prev = (u(2:end,:)+u(1:end-1,:))/2;     
%  plot((x_ar(2:end)),abs(u(1,:)./U_in_2))
%  hold all
       end
      u = u./U_in_2; 
      u_theor_svet(cur_f_n) = u(1,end);
      u_theor_tenj(cur_f_n) = u(64,end);
end
        
%  %% plotting     
 figure;
   plot(f_ar,abs(u_theor_svet))
  hold all
  plot(f_ar,abs(u_theor_tenj))
%  plot((x_ar(2:end)),abs(u(1,:)))
%  hold all
%  plot((x_ar(2:end)),abs(u(64,:)))

save('theor_field_01','u_theor_svet','u_theor_tenj','f_ar','rho_in','x_in')
