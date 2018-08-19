% hold on;
% close all;
figure;
Phi = 0/180*pi;

T0 =   0/180*pi;


T1 =  45/180*pi;


number = 100;

curve_number = 5;




symcase = 'asym';

%%

x = linspace(0,2,number);

temp1 = cos(T0+(T1-T0).*x);

temp2 = -log(temp1) + log(cos(T0));

y = 1/(T1-T0).*temp2;

plot(x,y,'r-');hold on;
plot(-x,-y,'r-');hold on;





switch symcase
    
    case 'asym'
        
        plot(-x,y,'b-');hold on;
        plot(x,-y,'b-');hold on;
end

% % cord =[ x',y'];
% % 
% % [x_rotated,y_rotated] = Rz_group(cord,Phi,[0 0 0]);
% % 
% %  plot(x_rotated,y_rotated,'b-');hold on;
% %  plot(-x_rotated,-y_rotated,'b-');hold on;
% % 
% % 
% % switch symcase
% %     
% %     case 'asym'
% %         
% %         plot(-x_rotated,y_rotated,'g-');hold on;
% %         plot(x_rotated,-y_rotated,'g-');hold on;
% % end
% % 


%%

for ii = 1:curve_number
    
    plot(x,y+1/(curve_number/2)*ii,'r-');hold on;
    plot(x,y-1/(curve_number/2)*ii,'r-');hold on;
    
    
    plot(-x,-(y+1/(curve_number/2)*ii),'r-');hold on;
    plot(-x,-(y-1/(curve_number/2)*ii),'r-');hold on;
    
    
    
    switch symcase
        
        case 'asym'
            
            
            
            plot(-x,y+1/(curve_number/2)*ii,'b-');hold on;
            plot(-x,y-1/(curve_number/2)*ii,'b-');hold on;
            
            
            plot(x,-(y+1/(curve_number/2)*ii),'b-');hold on;
            plot(x,-(y-1/(curve_number/2)*ii),'b-');hold on;
    end
    
    
    
    
end




axis image;axis([-1,1,-1,1]);

xlabel('Normalized length direction, x')
ylabel('Normalized width direction')
grid on;

