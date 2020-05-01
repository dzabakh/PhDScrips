clear all

figure;
hold all ;
N_x = 70;
N_y = 70;
N_z = 70;

c = 343;

n_t = 300;
dx = 0.02;
dt = 1/40000;


num_nodes = N_x*N_y*N_z;


%% calculating cordinates

node_coords = zeros(num_nodes , 3) ;
cur_num_node = 0 ; 
for z_idx = 0:N_z -1
   for y_idx = 0: N_y - 1 
       for x_idx = 0:N_x - 1
            cur_num_node = cur_num_node+1 ;  
            node_coords(cur_num_node , :) = dx*[x_idx,y_idx,z_idx];
        end
   end 
end
%% Mach_matrix

  M = zeros(num_nodes,1);
  R_node = 0*M;
  
  for cur_node=1: num_nodes
    R_node(cur_node) = sqrt((node_coords(cur_node,1)-0.7).^2 + ...
        (node_coords(cur_node,3)-0.7).^2); 
  end
  V_array = [0,20,40,60,80] ;
  
 for v_ind = 1 
     V = V_array(v_ind) ;
%      m = V/c ;
%      M =  m*0.5*(1 - tanh((R_node - 0.2)/0.02)) ;
%      M = Mach(node_coords(:,1) - dx*N_x/2, node_coords(:,2), ...
%          node_coords(:,3)-dx*N_z/2, V, 0.2, dx) ;
%      
% 
% %     assembling of the K matrix
%     K = sparse(num_nodes,num_nodes);
% 
%     %first edge of rectangle
%     K(1,1) = -6;
%     K(1,2) = 1;
%     K(1,1+N_x) = 1;
% 
%     %second edge of rectangle
%     K(N_x,N_x) = -6;
%     K(N_x,N_x-1) = 1;
%     K(N_x,N_x+N_x) = 1;
% 
%     %third edge of rectangle
%     K(end-N_x+1,end-N_x+1) = -6;
%     K(end-N_x+1,end-N_x+2) = 1;
%     K(end-N_x+1,end-2*N_x+1) = 1;
% 
%     %fourth edge of rectangle
% 
%     K(end,end) = -6;
%     K(end,end-1) = 1;
%     K(end,end-N_x) = 1;
% 
%     p0 = waitbar(0,'assembling K matrix');
% 
%     for cur_num = 1: num_nodes
%          waitbar(cur_num/num_nodes)
%           n_x = node_coords(cur_num,1)/dx+1;
%           n_y = node_coords(cur_num,2)/dx+1;
%           n_z = node_coords(cur_num,3)/dx+1;
% 
%           n_xnext = cur_num + 1;
%           n_ynext = cur_num + N_x;
%           n_znext = cur_num + N_y*N_x;
% 
%           n_xprev = cur_num - 1;
%           n_yprev = cur_num - N_x;
%           n_zprev = cur_num - N_y*N_x;
% 
% 
%         if(n_x>1)&&(n_y>1)&&(n_z>1)&&(n_x<N_x)&&(n_y<N_y)&&(n_z<N_z)
%             K(cur_num,cur_num) = -6+2*(M(cur_num))^2;
%             K(cur_num,n_xprev)= 1;
%             K(cur_num,n_xnext)= 1;
%             K(cur_num,n_yprev)= 1-(M(n_yprev))^2;
%             K(cur_num,n_ynext)= 1-(M(n_ynext))^2;
%             K(cur_num,n_zprev)= 1;
%             K(cur_num,n_znext)= 1;
%         elseif (n_x==1)&&(n_y>1)&&(n_z>1)&&(n_x<N_x)&&(n_y<N_y)&&(n_z<N_z) 
%             K(cur_num,cur_num) = -6+2*(M(cur_num))^2;
%             K(cur_num,n_xnext)= 1;
%             K(cur_num,n_yprev)= 1-(M(n_yprev))^2;
%             K(cur_num,n_ynext)= 1-(M(n_ynext))^2;
%             K(cur_num,n_zprev)= 1;
%             K(cur_num,n_znext)= 1;
%         elseif (n_x>1)&&(n_y==1)&&(n_z>1)&&(n_x<N_x)&&(n_y<N_y)&&(n_z<N_z) 
%             K(cur_num,cur_num) = -6+2*(M(cur_num))^2;
%             K(cur_num,n_xprev)= 1;
%             K(cur_num,n_xnext)= 1;
%             K(cur_num,n_ynext)= 1-(M(n_ynext))^2;
%             K(cur_num,n_zprev)= 1;
%             K(cur_num,n_znext)= 1;
%         elseif (n_x>1)&&(n_y>1)&&(n_z==1)&&(n_x<N_x)&&(n_y<N_y)&&(n_z<N_z) 
%             K(cur_num,cur_num) = -6+2*(M(cur_num))^2;
%             K(cur_num,n_xprev)= 1;
%             K(cur_num,n_xnext)= 1;
%             K(cur_num,n_yprev)= 1-(M(n_yprev))^2;
%             K(cur_num,n_ynext)= 1-(M(n_ynext))^2;
%             K(cur_num,n_znext)= 1;
%         elseif (n_x>1)&&(n_y>1)&&(n_z>1)&&(n_x==N_x)&&(n_y<N_y)&&(n_z<N_z) 
%             K(cur_num,cur_num) = -6+2*(M(cur_num))^2;
%             K(cur_num,n_xprev)= 1;
%             K(cur_num,n_yprev)= 1-(M(n_yprev))^2;
%             K(cur_num,n_ynext)= 1-(M(n_ynext))^2;
%             K(cur_num,n_zprev)= 1;
%             K(cur_num,n_znext)= 1;
%         elseif (n_x>1)&&(n_y>1)&&(n_z>1)&&(n_x<N_x)&&(n_y==N_y)&&(n_z<N_z) 
%             K(cur_num,cur_num) = -6+2*(M(cur_num))^2;
%             K(cur_num,n_xprev)= 1;
%             K(cur_num,n_xnext)= 1;
%             K(cur_num,n_yprev)= 1-(M(n_yprev))^2;
%             K(cur_num,n_zprev)= 1;
%             K(cur_num,n_znext)= 1;
%         elseif (n_x>1)&&(n_y>1)&&(n_z>1)&&(n_x<N_x)&&(n_y<N_y)&&(n_z==N_z) 
%             K(cur_num,cur_num) = -6+2*(M(cur_num))^2;
%             K(cur_num,n_xprev)= 1;
%             K(cur_num,n_xnext)= 1;
%             K(cur_num,n_yprev)= 1-(M(n_yprev))^2;
%             K(cur_num,n_ynext)= 1-(M(n_ynext))^2;
%             K(cur_num,n_zprev)= 1;
%         end
%     end
%     close (p0)
% 
%     % assembling of Dy matrix
%     Dy = sparse(N_x*N_y,N_x*N_y);
% 
%     for cur_num= 1: num_nodes
% 
%         n_y = node_coords(cur_num,2)/dx+1; %% coordinates on y-axis
% 
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
% 
% 
%       K = c^2*K/dx^2;
%       Dy = 2*c*Dy/(2*dx);
% 
%      save(['real_KD_matricies_m70',num2str(V)],'K','Dy')
      load(['KD_matricies_m70',num2str(V)],'K','Dy')



    u_cur = zeros(num_nodes,1);
    v_cur = zeros(num_nodes,1);

     u_cur(1:num_nodes,1) = 0*exp( - ((node_coords(:,1)-0.3).^2 + (node_coords(:,2)-0.5).^2+(node_coords(:,3)-0.5).^2) *155);

    for n= 1:size(node_coords,1)
       if( abs(node_coords(n,1) - 0.4) < dx ) && ( abs(node_coords(n,2) - 1) < dx) && (abs(node_coords(n,3) - 0.7) < dx ) 
            source_idx = n   
       end
    end    
    if source_idx == 0 
        'trouble with source'
    end

    for n= 1:size(node_coords,1)
       if( abs(node_coords(n,1) - 1) < dx ) && ( abs(node_coords(n,2) - 0.4) < dx ) && (abs(node_coords(n,3) - 0.7) < dx ) 
            reciver_idx_1 = n   
       end
    end

    for n= 1:size(node_coords,1)
       if( abs(node_coords(n,1) - 1) < dx ) && ( abs(node_coords(n,2) - 0.7) < dx )&& (abs(node_coords(n,3) - 0.7) < dx )  
            reciver_idx_2 = n   
       end
    end

     idx_source_arr = zeros(num_nodes,1);
     idx_source_arr(source_idx) =1;

    idx_reciver_arr = zeros(num_nodes,1);
    idx_reciver_arr(reciver_idx_1) =1;


    source_t = zeros(1,n_t);
    source_t(1) = 1;
    %% Ruunge kutta
    res = zeros(num_nodes,n_t);

     p = waitbar(0,'bambuk');

    for num_t = 1:n_t
        waitbar(num_t/n_t)
        k1u = v_cur; 
        k1v= K*u_cur - Dy*v_cur ;

        k2u = v_cur+ dt/2*k1v;
        k2v = K*(u_cur+dt/2*k1u)- Dy*(v_cur+dt/2*k1v);

        k3u = v_cur+ dt/2*k2v;
        k3v = K*(u_cur+dt/2*k2u)- Dy*(v_cur+dt/2*k2v);

        k4u = v_cur+ dt*k3v;
        k4v = K*(u_cur+dt*k3u)- Dy*(v_cur+dt*k3v);

        v_cur = v_cur + dt/6*(k1v+2*k2v+2*k3v+k4v) + dt*source_t(num_t)*idx_source_arr; 
        u_cur = u_cur + dt/6*(k1u+2*k2u+2*k3u+k4u);

    %      max(u_cur)

         res(:,num_t) = u_cur;
    end
    close(p)

    %%
    
    
    res_p = zeros(N_x*(N_y-1)*N_z, n_t);
    for i = 1:(N_y*N_z - 1)
        res_p((N_x*(i-1) + 1):N_x*i,:) = ...
            (res((N_x*i + 1):N_x*(i+1),:) - res((N_x*(i-1) + 1):N_x*i,:))/dx ;
    end
    
    
    res_p_r1 = res_p(reciver_idx_1,:);
    res_r1 = res(reciver_idx_1,:);
    res_r2 = 0;%res(reciver_idx_2,:);
    save(['res70', num2str(V), '.mat'], 'res_r1', 'res_r2')
    
    res_r1_f = ifft(res_r1);
    res_r2_f = ifft(res_r2);
    res_p_r1_f = ifft(res_p_r1) ;
    
    sz = length(res_r1);
    f_max = 1/dt;
    f_array = (0:sz-1)/(dt*sz) ;


    f_high = 3000;
    f_width = 200;

    mask_1 = 0.5*(1-tanh((f_array - f_high)/f_width)) + 0.5*(1+tanh((-f_max + f_array + f_high)/f_width)) ; 

    f_low = 200 ; 
    f_width = 200 ; 

    mask_2 = (0.5*(1+tanh((f_array - f_low)/f_width))) .*( 0.5*(1-tanh((-f_max + f_array + f_low)/f_width))) ; 

    mask = mask_1 .* mask_2 ;


    res_p_r1_filt = real(fft(res_p_r1_f.*mask));
    res_r1_filt = real(fft(res_r1_f.*mask));
    res_r2_filt = real(fft(res_r2_f.*mask));


    plot(c*(1:n_t)*dt, -res_p_r1_filt/1.506e-11)
%     plot(c*(1:n_t)*dt, res_r2_filt)
    % save('res3d','res','N_x','N_y')

 end
