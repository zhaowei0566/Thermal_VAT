% hold on;
close all;
figure;
Phi = 0/180*pi;

T0 =  60 /180*pi;


T1 =  0/180*pi;




number =24;
%%
x0 = 0;

Y0 = linspace(0,1,number );


symcase = 'sym';

for ii = 1:length(Y0)
   
    
    
    y0 = Y0(ii);
    
    Theta=Phi+(T1-T0)*x0+T0;
    
    k = tan(Theta);
    
    c = y0-k*x0;
    
    x = linspace(0,1,100);
    
    y = tan(Phi+(T1-T0).*x+T0).*x + c;
    
    plot(x,y,'r-');hold on;
    plot(-x,-y,'r-');hold on;
    
    
    switch symcase
        
        case 'asym'
            plot(-x,y,'k-');hold on;
            plot(x,-y,'k-');hold on;
    end
end

axis([-1 1 -1 1])

%%
y0 = -1;

X0 = linspace(0,1,number /2);



for ii = 1:length(X0)
   
    
    
    x0 = X0(ii);
    
    Theta=Phi+(T1-T0)*x0+T0;
    
    k = tan(Theta);
    
    c = y0-k*x0;
    
    x = linspace(0,1,100);
    
    y = tan(Phi+(T1-T0).*x+T0).*x + c;
    
    plot(x,y,'r-');hold on;
    plot(-x,-y,'r-');hold on;
    
    
    switch symcase
        
        case 'asym'
            plot(-x,y,'k-');hold on;
            plot(x,-y,'k-');hold on;
    end

    
end

%%


