clear ;
load('35mps');
c0 = 343 ;
U_m = [U_m_0, U_m_300, U_m_600].' ;

Z_real = [890; 590; 290] ;

U_forplot = rot90(rot90(U_m_0(1:floor(end/2)))) ;

figure;
plot(Y(1:floor(end/2)), U_forplot/c0);

%%
% M_600 = max(U_m_600)/c0 * 0.48 .* (1 - tanh((Y - 290)/20)) ;
% M_300 = max(U_m_300)/c0 * 0.48 .* (1 - tanh((Y - 290)/35)) ;
% M_0 = max(U_m_0)/c0 * 0.48 .* (1 - tanh((Y - 290)/50)) ;

hold all;
plot(Y, M_600) ;

% Podobral parametry approxim. func. dlya profilya 
% skorosti v kazhdom sluchae. Shirina yadra perenositsya s 
% koefficientom 1/1.5, shirina pogran. sloya s 1/\sqrt(1.5), 
% poskol'ku zavisit ot skorosti potoka i \nu.
% Vse beretsya is podobiya po Re.

%%
% x = 0:0.05:0.6 ;
% p = polyfit([0.290, 0.590, 0.890], [20, 35, 50], 2) ;
% p2 = polyval(p, x) ;
% figure; plot(x,p2) ;

% Shirina yadra spadaet lineyno - OK
%%