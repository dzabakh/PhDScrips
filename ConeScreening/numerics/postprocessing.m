load('tenj')

sp_tenj = sig_int_sp;

r_tenj = r;

load('svet')

sp_svet = sig_int_sp;

r_svet = r;

dr1 = 0.008;%0.01;
dr2 = 0.0025;%0.015;

load('theor_field_01')

n_exp = 1;

k = f_ar*2*pi/343;

%%
figure;
plot(f_ar,real(1i*log(sp_svet(n_exp,188:62:625)./u_theor_svet)./k))
 hold all
 plot(f_ar,real(1i*log(sp_tenj(n_exp,188:62:625)./u_theor_tenj)./k))
% hold all
% plot(f_ar,real(u_theor_svet),'b')
% hold all
% plot(f_ar,real(u_theor_tenj),'k')
% 
% figure;
% plot(f_array,abs(sp_svet(n_exp,:).*exp(1i*f_array*2*pi/343*dr)),'b.')
% hold all
% plot(f_array,abs(sp_tenj(n_exp,:).*exp(1i*f_array*2*pi/343*dr)),'k.')
% hold all
% plot(f_ar,abs(u_theor_svet),'b')
% hold all
% plot(f_ar,abs(u_theor_tenj),'k')

%%
figure;
plot(f_array,real(sp_svet(n_exp,:).*exp(1i*f_array*2*pi/343*dr1)))
hold all
plot(f_array,real(sp_tenj(n_exp,:).*exp(1i*f_array*2*pi/343*dr2)))
hold all
plot(f_ar,real(u_theor_svet),'b')
hold all
plot(f_ar,real(u_theor_tenj),'k')

figure;
plot(f_array,imag(sp_svet(n_exp,:).*exp(1i*f_array*2*pi/343*dr1)))
hold all
plot(f_array,imag(sp_tenj(n_exp,:).*exp(1i*f_array*2*pi/343*dr2)))
hold all
plot(f_ar,imag(u_theor_svet),'b')
hold all
plot(f_ar,imag(u_theor_tenj),'k')
