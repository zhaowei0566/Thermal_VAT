% hold on;
% close all;
figure;

t0 =45;
t1 = -45;
T0T1 = [ t0 t1 ;
    -t0 -t1;
    t0 t1;
    -t0 -t1
    -t0 -t1;
    t0 t1;
    -t0 -t1;
    t0 t1];


colorss=['krbggbrk'];

for ee = 1:size(T0T1,1)
    
    
    Phi = 0/180*pi;
    
    T0 =  T0T1(ee,1)/180*pi;
    
    
    T1 =  T0T1(ee,2)/180*pi;
    
    
    number = 100;
    
    curve_number = 20;
    
    
    z = ee/size(T0T1,1)*ones(number,1);
    
    %%
    
    x = linspace(0,2,number);
    
    temp1 = cos(T0+(T1-T0).*x);
    
    temp2 = -log(temp1) + log(cos(T0));
    
    y = 1/(T1-T0).*temp2;
    
    plot3(x,y,z,[colorss(ee) '-'],'LineWidth',1);hold on;
    plot3(-x,-y,z,[colorss(ee) '-'],'LineWidth',1);hold on;
    
    
    %%
    
    for ii = 1:curve_number
        
        plot3(x,y+1/(curve_number/2)*ii,z,[colorss(ee) '-'],'LineWidth',1);hold on;
        plot3(x,y-1/(curve_number/2)*ii,z,[colorss(ee) '-'],'LineWidth',1);hold on;
        
        
        plot3(-x,-(y+1/(curve_number/2)*ii),z,[colorss(ee) '-'],'LineWidth',1);hold on;
        plot3(-x,-(y-1/(curve_number/2)*ii),z,[colorss(ee) '-'],'LineWidth',1);hold on;
        
        
    end
    
    
    
end


axis([-1,1,-1,1]);

xlabel('Normalized length direction, x')
ylabel('Normalized width direction')
grid on;

