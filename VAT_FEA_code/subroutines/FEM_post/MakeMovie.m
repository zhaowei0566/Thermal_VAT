%%
%% Test
% Test
% clear; close all;
%
% % Creates a 2D Mesh to plot surface
% x=linspace(0,1,100);
% [X,Y] = meshgrid(x,x);
% DELAY = 0.05;
% N=100; % Number of frames
% for i = 1:N
%     % Example of plot
%     Z = sin(2*pi*(X-i/N)).*sin(2*pi*(Y-i/N));
%
%     subplot(1,2,1)
%     surf(X,Y,Z);view(2);
%
%     subplot(1,2,2)
%
%     plot(x(i),x(i)^2,'ro');axis([0 1 0 1]);hold on;
%
%     pause(DELAY);
%
%     % Store the frame
%     M(i)=getframe(gcf); % leaving gcf out crops the frame in the movie.
% end
%
% % Output the movie as an avi file
% movie2avi(M,'Surface_movie.avi');

%%


% x1 = linspace(0,1);
% x2 = linspace(3/4,1);
% y1 = sin(2*pi*x1);
% y2 = sin(2*pi*x2);
%
% figure(2)
%
% % plot on large axes
% plot(x1,y1)
%
% % create smaller axes in top right, and plot on it
% axes('Position',[.7 .7 .2 .2])
% box on
% plot(x2,y2)



%% Make move in mode shape change with in-plane load factor
close all;!del *.avi
DELAY = 0.1;
Lambda_b=[-1:0.005:1];%0.7 0.8:0.01:1]%[-1:0.01:-0.8 -0.7:0.1:1]
for jj=1:length(Lambda_b)
    %   for lambda_b_ratio=[0:0.005:0.1 0.2:0.1:0.7 0.8:0.01:1]
    %     lambda_b_ratio=0.8;
    %
    lambda_b_ratio=Lambda_b(jj);
    Keff=Ktotal+lambda_b_ratio*buckLF*KGtotal; %% elastic stiffness + geometric stiffness
    
    [V,D]=eigs(Keff(ActiveDof,ActiveDof),Mtotal(ActiveDof,ActiveDof),10,'sm');
    
    %%
    Keff=Ktotal+lambda_b_ratio*buckLF*KGtotal; %% ela
    %%
    
    
    [DD,modeNo]=sort(diag(D));
    
    loadfactor=DD(1:10);
    cycleFreq=sqrt(loadfactor);
    VVsort=V(:,modeNo);
    frequency=cycleFreq/2/pi;
    vibrparameter(jj,1:4)=cycleFreq(1:4)*Stru.length^2/Stru.thickness*sqrt(Mat.density/Mat.E2);
    
    
    modeshape4movie;
    
    
    % only for in-plane bending load due to the numerical error
    if jj==length(Lambda_b)
        
        vibrparameter(jj,1)=0;
    end
    
    subplot(1,2,2);hold on;
    plot(lambda_b_ratio,vibrparameter(jj,1)/18.9,'MarkerFaceColor',[1 0 0],...
        'MarkerEdgeColor',[1 0 0],...
        'MarkerSize',8,...
        'Marker','o',...
        'LineStyle','none',...
        'Color',[1 0 0]);
    axis([-1 1 0 1]);box on;grid on;
    xlabel('In-plane load factor','FontSize',12);
    ylabel('Fequency ratio, \lambda_1/\lambda_n','FontSize',12);
    
    
    pause(DELAY);
    M(jj)=getframe(gcf);
end


movie2avi(M,'ShearPreVibration.avi','fps',15);
