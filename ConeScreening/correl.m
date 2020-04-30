function [corrf] = correl(sig1, sig2)

    m = min([length(sig2), length(sig1)]) ;
    sig2 = sig2(1:m) ;
    sig1 = sig1(1:m) ;
    corrf = (fft(ifft(sig1) .* conj(ifft(sig2)))) ;
end    
    

