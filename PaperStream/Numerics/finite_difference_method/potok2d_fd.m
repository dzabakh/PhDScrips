clear all
c = 343/2 ;

N_x = 200;
N_y = 200;
r_nozzle = 0.2 ;
m = 0/c; %% mach number
rho = 1.21 ;

r_source = [0, 0] ;
r_mic = [0, 0.3] ;

n_t = 300 ;
dt = 1/48000;
dx = 0.02;

num_nodes = N_x*N_y;

x_source = floor(r_source(1)/dx) + N_x/2 ;
y_source = floor(r_source(2)/dx) + N_y/2 ;

x_mic = floor(r_mic(1)/dx) + N_x/2 ;
y_mic = floor(r_mic(2)/dx) + N_y/2 ;

%% point source
t = [dt:dt:n_t*dt] ;
% y = gaussmf(t, [3*dt 15*dt]) ;

t_high = dt*15 ; 
t_width =  dt*4 ;
y = 1 - 0.5*(1-tanh((t - t_high)/t_width));

ys = (y(2:end) - y(1:end-1))/dt ;
ys = [ys, 0] ;

%% calculating cordinates
node_coords = zeros(num_nodes , 2) ;
    cur_num_node = 0 ; 

   for y_idx = 0: N_y - 1 
       for x_idx = 0:N_x - 1
            cur_num_node = cur_num_node + 1 ;  
            node_coords(cur_num_node , :) = [x_idx, y_idx] ;
        end
   end 

%% Mach_matrix

%  M = zeros(num_nodes,1);

%  up_coord   = node_coords(:,1) > 30 ;
%  down_coord = node_coords(:,1) < 70 ;
%  radial_m = sqrt((node_coords(:,1)-N_x/2).^2 + ...
%      (node_coords(:,2)-N_y/2).^2) < r_nozzle/dx ;
 
M = ones(num_nodes, 1) ;

Mxy = m * ...
(1 - tanh((sqrt((node_coords(:,1)-N_x/2).^2 + (node_coords(:,2)-N_y/2).^2)...
- r_nozzle/dx)/r_nozzle));


M = 0*M .* Mxy ;

% Me = reshape(M, N_x, N_y) ;
% figure; 
% imagesc(Me) ;

%  mask_stream = radial_m; %up_coord & down_coord ; 

%  M(mask_stream) = m ;



%% assembling of the K matrix

K = sparse(N_x * N_y, N_x * N_y) ;

% first edge of rectangle
K(1, 1)     = -4 ;
K(1, 2)     =  1 ;
K(1, 1+N_x) =  1 ;

%second edge of rectangle
K(N_x, N_x)     = -4 ;
K(N_x, N_x-1)   =  1 ;
K(N_x, N_x+N_x) =  1 ;

% third edge of rectangle
K(end-N_x+1,end-N_x+1)= -4;
K(end-N_x+1,end-N_x+2) = 1;
K(end-N_x+1,end-2*N_x+1) = 1;

% fourth edge of rectangle

K(end,end)=-4;
K(end,end-1) = 1;
K(end,end-N_x) = 1;

p = waitbar(0,'K matrix');
 
for cur_num= 1: N_x*N_y
    
    waitbar(cur_num/(N_x*N_y))
    
    n_j = ceil(cur_num/N_x); %% coordinates on y-axis
    
    n_i = cur_num - (n_j-1)*N_x; %% coordinates on x-axis
    
    n_iprev = cur_num-1;
    n_inext = cur_num+1;
    n_jprev = cur_num-N_x;
    n_jnext = cur_num+N_x;
    
    if(n_i>1)&&(n_j>1)&&(n_i<N_x)&&(n_j<N_y)
        K(cur_num,cur_num) = -4-(M(cur_num))^2;
        K(cur_num,n_iprev)= 1;
        K(cur_num,n_inext)= 1;
        K(cur_num,n_jprev)= 1-(M(n_jprev))^2;
        K(cur_num,n_jnext)= 1-(M(n_jnext))^2;
    elseif (n_j==1)&&(n_i>1)&&(n_i<N_x)
        K(cur_num,cur_num) = -4-(M(cur_num))^2;
        K(cur_num,n_iprev)= 1;
        K(cur_num,n_inext)= 1;
        K(cur_num,n_jnext)= 1-(M(n_jnext))^2;
    elseif (n_j==N_y)&&(n_i>1)&&(n_i<N_x)   
        K(cur_num,cur_num) = -4-(M(cur_num))^2;
        K(cur_num,n_iprev)=1;
        K(cur_num,n_jprev)= 1-(M(n_jprev))^2;
        K(cur_num,n_inext)= 1;
    elseif (n_i==1)&&(n_j>1)&&(n_j<N_y)
        K(cur_num,cur_num) = -4-(M(cur_num))^2;
        K(cur_num,n_inext)= 1;
        K(cur_num,n_jprev)= 1-(M(n_jprev))^2;
        K(cur_num,n_jnext)= 1-(M(n_jnext))^2;
    elseif (n_i==N_x)&&(n_j>1)&&(n_j<N_y)
        K(cur_num,cur_num) = -4-(M(cur_num))^2;
        K(cur_num,n_iprev)=1;
        K(cur_num,n_jprev)= 1-(M(n_jprev))^2;
        K(cur_num,n_jnext)= 1-(M(n_jnext))^2;
    end
end

close(p) ;
K = c^2*K/dx^2;

%% assembling of Dy matrix
Dy = sparse(N_x*N_y,N_x*N_y);


for cur_num= 1: N_x*N_y

    n_j = ceil(cur_num/N_x); %% coordinates on y-axis
    
    
    n_jprev = cur_num-N_x;
    n_jnext = cur_num+N_x;
    
    if (n_j>1)&&(n_j<N_y)
       Dy(cur_num,n_jnext) = 1*M(n_jnext);
       Dy(cur_num,n_jprev) = -1*M(n_jprev);
    elseif (n_j==1)
       Dy(cur_num,n_jnext) = 1*M(n_jnext); 
    elseif (n_j==N_y)
       Dy(cur_num,n_jprev) = -1*M(n_jprev); 
    end
end

Dy = 2*c*Dy/(2*dx);

%% Ruunge kutta
u_cur = zeros(num_nodes,1);
v_cur = zeros(num_nodes,1);

% u_cur = 1/10*exp( - ((node_coords(:,1)-50).^2 + (node_coords(:,2)-50).^2) / 2^2 );
% u_cur(x_source + y_source*N_y, 1) = 1 ;

res = zeros(N_x*N_y,n_t);

p = waitbar(0,'Runge Kutta');


for num_t = 1:n_t
    waitbar(num_t/n_t)
    u_cur(x_source + y_source*N_y) = y(num_t) ;      
    k1u = v_cur; 
    k1v= K*u_cur - Dy*v_cur;
    
    k2u = v_cur+ dt/2*k1v;
    k2v = K*(u_cur+dt/2*k1u) - Dy*(v_cur+dt/2*k1v);
    
    k3u = v_cur+ dt/2*k2v;
    k3v = K*(u_cur+dt/2*k2u)- Dy*(v_cur+dt/2*k2v);
    
    k4u = v_cur+ dt*k3v;
    k4v = K*(u_cur+dt*k3u)- Dy*(v_cur+dt*k3v);
    
    v_cur = v_cur + dt/6*(k1v+2*k2v+2*k3v+k4v);
    u_cur = u_cur + dt/6*(k1u+2*k2u+2*k3u+k4u);
    
    max(u_cur)
    res(:,num_t) = u_cur;
end
close(p)
res_mic = res(x_mic + y_mic*N_y, :) ;
% res_mic = res_mic*4*pi/rho ;
save('res2d','res','res_mic','N_x','N_y','y','ys', 'dx')
