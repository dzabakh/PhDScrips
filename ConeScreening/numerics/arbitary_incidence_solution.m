clear all
k = 1;
a = 1;
ka2 = 1;
x0_ar= 0.2:0.5:12;
N_x = 1000;
N_add = 100;
h_im = 0.01 ;
N_iter = 1;
dphi = 0.05;
phi = 0:dphi:(2*pi-dphi);
phi0_ar = phi;
th_in = 0;
phi_in = pi;


ones_phi = ones(length(phi),1);






Kernel_full = zeros(length(x0_ar),length(phi0_ar),length(phi)-1,N_x+N_add-2);

u = [];
u_inv = [];
wb = waitbar(0 , 'smoking bamboo') ;
for n_x_cur = 1:length(x0_ar)
    
   for n_phi_cur = 1:length(phi0_ar)
    
        waitbar((n_x_cur-1)/(length(x0_ar)-1))

        x0 = x0_ar(n_x_cur);
        phi0 = phi0_ar(n_phi_cur); 
        dx = x0/N_x;
        dx_im = h_im/N_add;
        x_ar_sing_1 = 0 + 1i*(dx_im:dx_im:h_im);
        x_ar_sing_2 = 1i*h_im+(0+dx:dx:x0);
        x_ar_sing_3 = x0 + 1i*(h_im-dx_im:-dx_im:dx_im);
        x = [x_ar_sing_1,x_ar_sing_2,x_ar_sing_3];
%     x = dx:dx:(x0 - dx);
%     x_add = x(end):dx/N_add:(x0-dx/N_add);
%     x = [x,x_add];
    x_diff_grid = ones_phi.*(x0-x);
    x_grid = ones_phi*x;
    phi_grid = (phi.'-phi0)*ones(1,length(x));
    phi_in_grid = (phi.'-phi_in)*ones(1,length(x));
    
    U_in_grid= exp(-1i*k*x_grid.*(th_in^2/2+th_in*a*cos(phi_in_grid)));
    
    
    U_in_r= exp(-1i*k*x0.*(th_in^2/2+th_in*a*cos(phi_in - phi0)));
    
    dx_ar = (x(2:end)-x(1:end-1));
    dphi_ar = phi(2:end)-phi(1:end-1);


    Kernel = 1i*x_grid.*ka2/(2*pi).*(1./x_diff_grid + (x_grid-x0*cos(phi_grid))./x_diff_grid.^2)...
        .*exp(0.5i*ka2.*(x0^2+x_grid.^2-2*x_grid.*x0.*cos(phi_grid))./x_diff_grid).*((dphi_ar.')*dx_ar);
     Kernel = (Kernel(2:end,:)+Kernel(1:end-1,:))/2;
     Kernel = (Kernel(:,2:end)+Kernel(:,1:end-1))/2;
     
     Kernel_full(n_x_cur,n_phi_cur,:,:) = Kernel;

%       figure;
%       surf(abs(Kernel(1:20:end,1:20:end)))

%        for n_cur_it = 1:N_iter
% 
%         u_cur = dphi_ar*Kernel*dx_ar + 2*U_in_r;
%         end
%        u = [u,u_cur];
       
   end
             
end
close(wb)


 

%% here we compare kernel with the exact symmetrical
% 
% 
%  Kernel_sym = dphi_ar*Kernel;
% 
% 
% Kernel_sym_true = 1i*ka2*x0*x./(x0-x).^2.*exp(0.5i*ka2.*(x0^2+x.^2)./(x0-x))...
%     .*(besselj(0,ka2.*x0*x./(x0-x))+1i*besselj(1,ka2.*x0*x./(x0-x)));
% 
% Kernel_sym_true = (Kernel_sym_true(2:end)+Kernel_sym_true(1:end-1))/2;
%  
% 
% 
%  figure;
%  plot(x(2:end),abs(Kernel_sym(1:end)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 figure;
 plot(abs(ka2)*x0_ar,abs(u))