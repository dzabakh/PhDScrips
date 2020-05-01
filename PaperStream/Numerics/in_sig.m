dt = 0.01 ;
n_t = 500 ;
t_max = n_t*dt ;
t = [dt:dt:t_max] ;


t_high = dt*15 ; 
t_width =  dt*3 ;

mask = 1 - 0.5*(1-tanh((t - t_high)/t_width));
plot(t, mask) ;



