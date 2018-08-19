close all;
a= 2;

b = 2;


figure(23); 

x_RHS = linspace(0,a/2,101);

number_of_curves = 20;



theta = [0 0.1]

T0T1 = [theta
       -theta]
T0 = T0T1(1,1)/180*pi;
T1 = T0T1(1,2)/180*pi;

y = a/2/(T1-T0).*(-log(cos(T0+2*(T1-T0).*x_RHS/a))+log(cos(T0)));

subplot(2,2,1);
plot(x_RHS,y,'r-');hold on;
plot(-x_RHS,-y,'r-');axis([-a/2,a/2,-a/2,a/2]);
title('T_0=T_1 = 0^o');grid on



theta = [89.9 89.91]

T0T1 = [theta
       -theta]
T0 = T0T1(1,1)/180*pi;
T1 = T0T1(1,2)/180*pi;

y = a/2/(T1-T0).*(-log(cos(T0+2*(T1-T0).*x_RHS/a))+log(cos(T0)));

subplot(2,2,2);
plot(x_RHS,y,'r-');hold on;
plot(-x_RHS,-y,'r-');axis([-a/2,a/2,-a/2,a/2]);grid on;
title('T_0=T_1 = 90^o');



theta = [0 60]

T0T1 = [theta
       -theta]
T0 = T0T1(1,1)/180*pi;
T1 = T0T1(1,2)/180*pi;

y = a/2/(T1-T0).*(-log(cos(T0+2*(T1-T0).*x_RHS/a))+log(cos(T0)));

subplot(2,2,3);
plot(x_RHS,y,'r-');hold on;
plot(-x_RHS,-y,'r-');axis([-a/2,a/2,-a/2,a/2])
title('T_0=0^o, T_1 = 60^o');grid on



theta = [-30 60]

T0T1 = [theta
       -theta]
T0 = T0T1(1,1)/180*pi;
T1 = T0T1(1,2)/180*pi;

y = a/2/(T1-T0).*(-log(cos(T0+2*(T1-T0).*x_RHS/a))+log(cos(T0)));

subplot(2,2,4);
plot(x_RHS,y,'r-');hold on;
plot(-x_RHS,-y,'r-');axis([-a/2,a/2,-a/2,a/2])
title('T_0=-30^o, T_1 = 60^o');grid on




% 
% T0 = T0T1(2,1)/180*pi;
% T1 = T0T1(2,2)/180*pi;
% 
% 
% y = a/2/(T1-T0).*(-log(cos(T0+2*(T1-T0).*x_RHS/a))+log(cos(T0)));
% 
% figure(23); 
% plot(x_RHS,y,'r-');hold on;
% plot(-x_RHS,-y,'r-');
% 



axis([-a/2,a/2,-a/2,a/2])




