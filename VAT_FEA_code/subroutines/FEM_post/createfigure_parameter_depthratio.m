function createfigure_parameter_depthratio(X1, YMatrix1)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 11-Oct-2015 22:39:04

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YMinorTick','on','XMinorTick','on',...
    'LineWidth',2,...
    'FontSize',16);
%% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[10 80]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1,'LineWidth',3);
set(plot1(1),'DisplayName','Design A, [(0/90)_2]_s');
set(plot1(2),'DisplayName','Design B, [(0/90)_2]_s','Color',[0 0 0]);
set(plot1(3),'Color',[1 0 0],'DisplayName','Design C, [(0/90)_2]_s');
set(plot1(4),'LineStyle','--','DisplayName','Design A, [(45/-45)_2]_s');
set(plot1(5),'LineStyle','--','DisplayName','Design B, [(45/-45)_2]_s',...
    'Color',[0 0 0]);
set(plot1(6),'LineStyle','--','Color',[1 0 0],...
    'DisplayName','Design C, [(45/-45)_2]_s');

% Create xlabel
xlabel('Stiffener depth ratio, (h_s/b_s)','FontSize',16);

% Create ylabel
ylabel('Buckling parameter','FontSize',16);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'LineWidth',2);

