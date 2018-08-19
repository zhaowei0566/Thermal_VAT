% hold on;
% close all;
figure;
Phi = 0/180*pi;

T0 = -89./180*pi;


T1 = 45/180*pi;


number = 101;

curve_number = 20;




symcase = 'sym';

%%

physical_width = 0.75;
physical_length = 0.2;


x = linspace(0,physical_length,number);

% x = unique([linspace(0,0.05,10)  linspace(0.05,0.1,2)  linspace(0.1,0.15,2) linspace(0.15,0.2,10)]);

temp1 = cos(T0+(T1-T0).*x/(physical_length));

temp2 = -log(temp1) + log(cos(T0));

y = 1/(T1-T0).*temp2;




plot(x,y,'r-');hold on;
plot(-x,-y,'r-');hold on;





switch symcase
    
    case 'asym'
        
        plot(-x,y,'k-');hold on;
        plot(x,-y,'k-');hold on;
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
            
            
            
            plot(-x,y+1/(curve_number/2)*ii,'k-');hold on;
            plot(-x,y-1/(curve_number/2)*ii,'k-');hold on;
            
            
            plot(x,-(y+1/(curve_number/2)*ii),'k-');hold on;
            plot(x,-(y-1/(curve_number/2)*ii),'k-');hold on;
    end
    
    
    
    
end




axis image;axis([-physical_length/2,physical_length/2,-physical_width/2,physical_width/2]);

xlabel('Normalized length direction, x')
ylabel('Normalized width direction')
grid on;

