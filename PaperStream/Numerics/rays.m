c = 343 ; 

L1 = 1.6 ;
L2 = 0.4 ; 

M=80 / 343 ; 

dz = -2.6 ; 

gg = (-0.9:0.001:0.9) ; 

gamma = gg / c ; 

fn = dz - L1 * (gamma ./ sqrt(1/c^2 - gamma.^2)) - L2*((gamma + (1/c- gamma*M)*M)./sqrt((1/c-gamma*M).^2-gamma.^2));

plot(gg, fn)

idx = 0 ;

for n = 1:length(gg)-1
   if (fn(n) >0 ) && (fn(n+1) <=0)
       idx = n ; 
   end
 
end

gamma_s = gg(idx) / c ; 

t = gamma_s* dz + sqrt(1/c^2 - gamma_s^2)* L1 + sqrt((1/c - gamma_s*M)^2 - gamma_s^2) * L2
t * c