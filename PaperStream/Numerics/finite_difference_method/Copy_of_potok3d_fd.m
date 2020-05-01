clear all

N_x = 70;
N_y = 70;
N_z = 70;

c = 343;
r_nozzle = 0.2 ;

r_source = [-0.3, 0.9, 0] ;
r_mic =    [0.3, 0.5, 0] ;

V_stream_ar = [0, 20, 40, 60, 80] ;



n_t = 80;
dt = 1/24000;
dx = 0.02;
rho0 = 1.21 ;
S = pi*(0.035/2)^2 ;

num_nodes = N_x * N_y * N_z ;

x_source = floor(r_source(1)/dx) + N_x/2 ;
y_source = floor(r_source(2)/dx) ;
z_source = floor(r_source(3)/dx) + N_z/2;

x_mic = floor(r_mic(1)/dx) + N_x/2 ;
y_mic = floor(r_mic(2)/dx) ;
z_mic = floor(r_mic(3)/dx) + N_z/2;

%% point source
t = [dt:dt:n_t*dt] ;
ys = t*0 ;
ys(3) = 1;
%% calculating cordinates

node_coords = zeros(num_nodes , 3) ;
    cur_num_node = 0 ; 
for z_idx = 0:N_z -1
  for y_idx = 0: N_y - 1 
      for x_idx = 0:N_x - 1
            cur_num_node = cur_num_node + 1 ;  
            node_coords(cur_num_node , :) = [x_idx,y_idx,z_idx];
       end
   end 
end
    
    
    
for v_ind = [1,3,5]
% %     
    V_stream = V_stream_ar(v_ind) ;    
    %% Mach_matrix
    M = Mach((node_coords(:,1) - N_x/2) * dx, ...
        (node_coords(:,2)) * dx, ...
        (node_coords(:,3)- N_z/2) * dx, V_stream, r_nozzle, dx) ;

    %% assembling of the K matrix

%     K = sparse(num_nodes,num_nodes);
%     % first edge of rectangle
%     K(1,1) = - 6;
%     K(1,2) =   1;
%     K(1,1+N_x) = 1;
%     %second edge of rectangle
%     K(N_x,N_x) = -6;
%     K(N_x,N_x-1) = 1;
%     K(N_x,N_x+N_x) = 1;
%     % third edge of rectangle
%     K(end-N_x+1,end-N_x+1)=-6;
%     K(end-N_x+1,end-N_x+2) = 1;
%     K(end-N_x+1,end-2*N_x+1) = 1;
%     % fourth edge of rectangle
%     K(end,end)=-6;
%     K(end,end-1) = 1;
%     K(end,end-N_x) = 1;
% 
%     p0 = waitbar(0,'Assembling K matrix');
% 
%     for cur_num= 1: num_nodes
%          waitbar(cur_num/num_nodes)
%           n_x = node_coords(cur_num,1)+1;
%           n_y = node_coords(cur_num,2)+1;
%           n_z = node_coords(cur_num,3)+1;
% 
%           n_xnext = cur_num+1;
%           n_ynext = cur_num + N_x;
%           n_znext = cur_num + N_y*N_x;
% 
%           n_xprev = cur_num-1;
%           n_yprev = cur_num - N_x;
%           n_zprev = cur_num - N_y*N_x;
% 
% 
%         if(n_x>1)&&(n_y>1)&&(n_z>1)&&(n_x<N_x)&&(n_y<N_y)&&(n_z<N_z)
%             K(cur_num,cur_num) = -6-(M(cur_num))^2;
%             K(cur_num,n_xprev)= 1;
%             K(cur_num,n_xnext)= 1;
%             K(cur_num,n_yprev)= 1-(M(n_yprev))^2;
%             K(cur_num,n_ynext)= 1-(M(n_ynext))^2;
%             K(cur_num,n_zprev)= 1;
%             K(cur_num,n_znext)= 1;
%         elseif (n_x==1)&&(n_y>1)&&(n_z>1)&&(n_x<N_x)&&(n_y<N_y)&&(n_z<N_z) 
%             K(cur_num,cur_num) = -6-(M(cur_num))^2;
%             K(cur_num,n_xnext)= 1;
%             K(cur_num,n_yprev)= 1-(M(n_yprev))^2;
%             K(cur_num,n_ynext)= 1-(M(n_ynext))^2;
%             K(cur_num,n_zprev)= 1;
%             K(cur_num,n_znext)= 1;
%         elseif (n_x>1)&&(n_y==1)&&(n_z>1)&&(n_x<N_x)&&(n_y<N_y)&&(n_z<N_z) 
%             K(cur_num,cur_num) = -6-(M(cur_num))^2;
%             K(cur_num,n_xprev)=1;
%             K(cur_num,n_xnext)= 1;
%             K(cur_num,n_ynext)= 1-(M(n_ynext))^2;
%             K(cur_num,n_zprev)= 1;
%             K(cur_num,n_znext)= 1;
%         elseif (n_x>1)&&(n_y>1)&&(n_z==1)&&(n_x<N_x)&&(n_y<N_y)&&(n_z<N_z) 
%             K(cur_num,cur_num) = -6-(M(cur_num))^2;
%             K(cur_num,n_xprev)= 1;
%             K(cur_num,n_xnext)= 1;
%             K(cur_num,n_yprev)= 1-(M(n_yprev))^2;
%             K(cur_num,n_ynext)= 1-(M(n_ynext))^2;
%             K(cur_num,n_znext)= 1;
%         elseif (n_x>1)&&(n_y>1)&&(n_z>1)&&(n_x==N_x)&&(n_y<N_y)&&(n_z<N_z) 
%             K(cur_num,cur_num) = -6-(M(cur_num))^2;
%             K(cur_num,n_xprev)=1;
%             K(cur_num,n_yprev)= 1-(M(n_yprev))^2;
%             K(cur_num,n_ynext)= 1-(M(n_ynext))^2;
%             K(cur_num,n_zprev)= 1;
%             K(cur_num,n_znext)= 1;
%         elseif (n_x>1)&&(n_y>1)&&(n_z>1)&&(n_x<N_x)&&(n_y==N_y)&&(n_z<N_z) 
%             K(cur_num,cur_num) = -6-(M(cur_num))^2;
%             K(cur_num,n_xprev)=1;
%             K(cur_num,n_xnext)= 1;
%             K(cur_num,n_yprev)= 1-(M(n_yprev))^2;
%             K(cur_num,n_zprev)= 1;
%             K(cur_num,n_znext)= 1;
%         elseif (n_x>1)&&(n_y>1)&&(n_z>1)&&(n_x<N_x)&&(n_y<N_y)&&(n_z==N_z) 
%             K(cur_num,cur_num) = -6-(M(cur_num))^2;
%             K(cur_num,n_xprev)=1;
%             K(cur_num,n_xnext)= 1;
%             K(cur_num,n_yprev)= 1-(M(n_yprev))^2;
%             K(cur_num,n_ynext)= 1-(M(n_ynext))^2;
%             K(cur_num,n_zprev)= 1;
%         end
%     end
%     close (p0)
%     K = c^2*K/dx^2;
%     save(['K_3d_dx002_70_V', num2str(V_stream)], 'K') ;
    load(['K_3d_dx002_70_V', num2str(V_stream)]) ;
% %     K = K * 0.05^2 / dx^2 ;
% %     assembling of Dy matrix
% 
%     Dy = sparse(N_x*N_y,N_x*N_y);
% 
%     for cur_num= 1: num_nodes
%         n_y = node_coords(cur_num,2)+1; %% coordinates on y-axis
%         n_yprev = cur_num-N_x;
%         n_ynext = cur_num+N_x;
% 
%         if (n_y>1)&&(n_y<N_y)
%            Dy(cur_num,n_ynext) = 1*M(n_ynext);
%            Dy(cur_num,n_yprev) = -1*M(n_yprev);
%         elseif (n_y==1)
%            Dy(cur_num,n_ynext) = 1*M(n_ynext); 
%         elseif (n_y==N_y)
%            Dy(cur_num,n_yprev) = -1*M(n_yprev); 
%         end
%     end
% 
%     Dy = 2*c*Dy/(2*dx);
%     save(['Dy2_dx002_70_3dV', num2str(V_stream)], 'Dy') ;
    
    load(['Dy2_dx002_70_3dV', num2str(V_stream)]) ; 
% %     Dy = Dy*0.05/dx ;
% % end
% 
% end
% % % Dy = Dy*0.01/dx ;
% % %% Ruunge kutta
u_cur = zeros(num_nodes,1);
v_cur = zeros(num_nodes,1);
u_sources = zeros(num_nodes,1);
% u_cur(1:num_nodes,1) = exp(-(...
%     (node_coords(:,1)-x_source).^2 + ...
%     (node_coords(:,2)-y_source).^2 + ...
%     (node_coords(:,3)-z_source).^2 ) / 8 );
% res = zeros(num_nodes,n_t);
p = waitbar(0,'Runge Kutta');
u_cur(x_source + y_source*N_x + z_source*N_x*N_y) = 1 ;
res_mic = t*0 ;

for num_t = 1:n_t
    waitbar(num_t/n_t)
    
    k1u = v_cur  ;
    k1v = K*u_cur - Dy*v_cur ;
    
    k2u = v_cur+ dt/2*k1v;
    k2v = K*(u_cur+dt/2*k1u) - Dy*(v_cur+dt/2*k1v) ;
    
    k3u = v_cur+ dt/2*k2v;
    k3v = K*(u_cur+dt/2*k2u) - Dy*(v_cur+dt/2*k2v);
     
    k4u = v_cur+ dt*k3v;
    k4v = K*(u_cur+dt*k3u)- Dy*(v_cur+dt*k3v) ;
    
    v_cur = v_cur + dt/6*(k1v+2*k2v+2*k3v+k4v);
    u_cur = u_cur + dt/6*(k1u+2*k2u+2*k3u+k4u);    
    
%     max(u_cur);
%     res(:,num_t) = u_cur;    
    
    res_mic(num_t) = u_cur(x_mic + y_mic*N_x + z_mic*N_x*N_y) ;    
end
close(p)
% res_mic = res(x_mic + y_mic*N_y + z_mic*N_z,:) ;
save(['res3d', num2str(V_stream)],'res_mic','N_x','N_y', 'ys')
end