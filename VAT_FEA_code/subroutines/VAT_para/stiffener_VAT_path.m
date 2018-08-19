close all;
a= 2;

b = 2;

x_RHS = linspace(0,a/2,100);

number_of_curves = 10;


T0 = T0T1(2,1)/180*pi;
T1 = T0T1(2,2)/180*pi;


y = a/2/(T1-T0).*(-log(cos(T0+2*(T1-T0).*x_RHS/a))+log(cos(T0)));



figure; hold on;


for curv = 1:number_of_curves
    
plot(x_RHS,y+(curv-1)*b/number_of_curves,'r-');hold on;

plot(x_RHS,y+(curv-1)*b/number_of_curves-a/2,'r-');

hold on;plot(-x_RHS,-(y+(curv-1)*b/number_of_curves),'r-');

hold on;plot(-x_RHS,-(y+(curv-1)*b/number_of_curves)+a/2,'r-');

end


axis([-a/2,a/2,-a/2,a/2])

%%

figure
% 
x_LHS = -fliplr(x_RHS);

y_LHS = a/2/(T1-T0).*(log(cos(T1+2*(T1-T0).*x_LHS/a))-log(cos(T1)));

% 
for curv = 1:number_of_curves
    hold on
    
plot(x_LHS,y_LHS+(curv-1)*b/number_of_curves); hold on;

plot(x_LHS,y_LHS+(curv-1)*b/number_of_curves-a/2);

hold on;plot(-x_LHS,-(y_LHS+(curv-1)*b/number_of_curves));

hold on;plot(-x_LHS,-(y_LHS+(curv-1)*b/number_of_curves)+a/2);

end


axis([-a/2,a/2,-a/2,a/2])

