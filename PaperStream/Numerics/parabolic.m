t_array = (-1:0.01:0.99) ; 

u = exp(-t_array.^2 /0.1^2) ;


A = 8000 ;

mask1 = exp( -A * (t_array + 1).^2  ) ;  
mask2 = exp( -A * (t_array - 1).^2  ) ;  

mask = mask1+ mask2 ; 

plot(mask)

sp = ifft(u) .* mask ; 

u_new = real(fft(sp)) ; 

plot(u_new)