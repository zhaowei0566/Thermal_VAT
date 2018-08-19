% %% Point move with a curve
% % DESCRIPTIVE TEXT
% %# control animation speed
% DELAY = 0.01;
% numPoints = 600;
% 
% %# create data
% x = linspace(0,10,numPoints);
% y = log(x);
% 
% %# plot graph
% figure('DoubleBuffer','on')                  %# no flickering
% plot(x,y, 'LineWidth',2), grid on
% xlabel('x'), ylabel('y'), title('y = log(x)')
% 
% %# create moving point + coords text
% hLine = line('XData',x(1), 'YData',y(1), 'Color','r', ...
%     'Marker','o', 'MarkerSize',6, 'LineWidth',2);
% hTxt = text(x(1), y(1), sprintf('(%.3f,%.3f)',x(1),y(1)), ...
%     'Color',[0.2 0.2 0.2], 'FontSize',8, ...
%     'HorizontalAlignment','left', 'VerticalAlignment','top');
% 
% %# infinite loop
% i = 1;                                       %# index
% while true      
%     %# update point & text
%     set(hLine, 'XData',x(i), 'YData',y(i))   
%     set(hTxt, 'Position',[x(i) y(i)], ...
%         'String',sprintf('(%.3f,%.3f)',[x(i) y(i)]))        
%     drawnow                                  %# force refresh
%     %#pause(DELAY)                           %# slow down animation
% 
%     i = rem(i+1,numPoints)+1;                %# circular increment
%     if ~ishandle(hLine), break; end          %# in case you close the figure
% end


%% SECTION TITLE
% DESCRIPTIVE TEXT
x = 1:100;
y = log(x);
DELAY = 0.05;
for i = 1:numel(x)
    clf;
    plot(x,y);
    hold on;
    plot(x(i),y(i),'r*');
    pause(DELAY);
    
    M(i)=getframe(gcf);
    
end

VideoWriter(M,'WaveMovie.avi');















