clear
clc



   primit = [5,17]; % 
%   primit = [5,8,10,13]; % 
%   primit=[3,13];%
% primit = [7,18];
% primit = [4,9,18,19];

M = primit(end);  % order of the sequence
L = 2^M -1; % length of sequence

num_primit = length(primit);

seed = zeros(1,M) ; 
seed(1) = 1;

result = zeros(L,1);

for n = 1:L
   result(n) = seed(1); 
   
   r = 0;
   for m = 1:num_primit
       r = xor(r,seed(primit(m)));
   end
   seed(2:M) = seed(1:(M-1));
   seed(1) = r ; 
      
end
result=0.5*(result-0.5);
result = [result;result;result];

ft=fft(result);
corr=real(ifft(ft.*conj(ft)));
figure;plot(corr);

SampleRate = 32768 ; 

filename = ['sample' num2str(M)];
 
wavwrite(result,SampleRate,16,filename)