function [x_arr] = my_sqrt(contour, x_0)
x_arr = contour * 0;
dx = 0.001;
    for n = 1:length(x_arr)
        for m = 1:5
            eta = x_0^2 - contour(n) ;
            mu = ((x_0 + dx)^2 - x_0^2) /dx ;
            x_0 = x_0 - eta /mu ;
        end
        x_arr(n) = x_0 ;
    end
end

