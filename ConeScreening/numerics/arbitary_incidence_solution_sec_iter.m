clear all

 
%% setting parameters

k = 1;
alpha = 1;
ka2 = k*alpha^2;
x_end = 1 ; 
N_x = 20;
n_im = 200;
h_im = 0.2;
n_pres = 1000;
N_iter = 20;
dx = x_end/N_x;
x_ar = dx:dx:x_end;
dphi = 0.1;
phi_ar = 0:dphi:(2*pi-dphi);
phi = 0;
th_in = 0;
phi_in = pi;


% %% calculating kernel
% 
% Kernel_full = zeros(length(phi_ar)-1,length(phi_ar),N_x,N_x);
% 
% 
% wb = waitbar(0 , 'smoking bamboo') ;
% for cur_n_1=1:length(x_ar)
%     x = x_ar(cur_n_1);
%     waitbar((cur_n_1-1)/(length(x_ar)-1))
%     for cur_n_phi = 1:length(phi_ar)
%         phi = phi_ar(cur_n_phi);
%         for cur_n_2 = 1:cur_n_1                      
%             if (cur_n_2==1)
%                 x3 = x_ar(cur_n_2+1);
%                 x2 = x_ar(cur_n_2);
%                 Kernel_full(:,cur_n_phi,cur_n_1,cur_n_2) = Kij_full(k,x,x3,x2,phi_ar,phi,alpha,n_pres,n_im,h_im,1);
%             elseif (cur_n_2 == cur_n_1)%&&(cur_n_2~=1)
%                 x2 = x_ar(cur_n_2);
%                 x1 = x_ar(cur_n_2-1);
%                 Kernel_full(:,cur_n_phi,cur_n_1,cur_n_2) = Kij_full(k,x,x2,x1,phi_ar,phi,alpha,n_pres,n_im,h_im,0);
%             else
%                 x3 = x_ar(cur_n_2+1);
%                 x2 = x_ar(cur_n_2);
%                 x1 = x_ar(cur_n_2-1);
%                 Kernel_full(:,cur_n_phi,cur_n_1,cur_n_2) = Kij_full(k,x,x2,x1,phi_ar,phi,alpha,n_pres,n_im,h_im,0)...
%                 + Kij_full(k,x,x3,x2,phi_ar,phi,alpha,n_pres,n_im,h_im,1);       
%             end
%         end
%     end
% end
% close(wb)
% Ker_add = Kernel_full(:,:,2:end,1);
% Kernel_full = Kernel_full(:,:,2:end,2:end);
% save('Kernel_full_1','Ker_add','Kernel_full')
%%  making iterations
load('Kernel_full_1')
% U_in_full = ones(length(phi_ar)-1,length(phi_ar)-1,N_x-1,N_x-1);
U_in_1 = ones(length(phi_ar)-1,N_x-1);
U_in_2 = ones(length(phi_ar),N_x-1);

u = 0*U_in_2;
u_prev = U_in_1;
figure;
       for n_cur_it = 1:N_iter
           for cur_n_x = 1:(length(x_ar)-1)
               for cur_n_phi = 1:(length(phi_ar)) 
                 u(cur_n_phi,cur_n_x) = sum(sum((squeeze(Kernel_full(:,cur_n_phi,cur_n_x,:)).*u_prev),1),2) + squeeze(sum(Ker_add(:,cur_n_phi,cur_n_x,:),1))+2.*U_in_2(cur_n_phi,cur_n_x);
               end
           end
           u_prev = (u(2:end,:)+u(1:end-1,:))/2;     
  plot(abs(ka2)*(x_ar(2:end)),abs(u(20,:)))
  hold all
        end
      
 

% %% plotting
% 
%  figure;
%  plot(abs(ka2)*(x_ar(2:end)),abs(u(20,:)))