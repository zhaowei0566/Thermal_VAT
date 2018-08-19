% get data out of MATLAB figure

D=get(gcf,'Children'); %get the handle of the line object
XData=get(D,'XData'); %get the x data
YData=get(D,'YData'); %get the y data
Data=[XData' YData']; %join the x and y data on one array nx2
%Data=[XData;YData]; %join the x and y data on one array 2xn


% % % % change the order of lines
chH = get(gcf,'Children');
set(gcf,'Children',[chH(end);chH(1:end-1)]);



% h = gcf; %current figure handle
% axesObjs = get(h, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
% xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
% ydata = get(dataObjs, 'YData');
% zdata = get(dataObjs, 'ZData');
