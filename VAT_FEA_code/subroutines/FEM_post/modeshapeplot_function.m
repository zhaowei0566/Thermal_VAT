
% addpath('V:\Wei_Zhao\OpenProgram\export_fig');

FEM.frequency=frequency;
FEM.modeshape=VVsort;

bendingmode=find(ActiveDof<=FEM.NodeNumber);

ActiveBendDOF=ActiveDof(bendingmode);

bendModeNo=length(bendingmode);

%---------------Mode Shape Plot-------------------

% if bendModeNo<10
%     mode=bendModeNo;
% else
%     mode=plotmodenNo;
% end

X=zeros(size(FEM.elementNodes,1),size(FEM.nodesCord,2));
Y=zeros(size(FEM.elementNodes,1),size(FEM.nodesCord,2));
Z=zeros(size(FEM.elementNodes,1),size(FEM.nodesCord,2));

coordinates=FEM.nodesCord;
nodes(1:size(FEM.elementNodes,1),1)=1:size(FEM.elementNodes,1);
nodes(1:FEM.elementNumber,2:size(FEM.elementNodes,2)+1)=FEM.elementNodes;
%         %---Plot Mesh----
% PlotMesh(coordinates,nodes)
% view(2)
%

deformUX=zeros(FEM.NodeNumber,1);
deformUY=zeros(FEM.NodeNumber,1);
deformUZ=zeros(FEM.NodeNumber,1);



for modeNo= real_mode_number
    
    deformUZ(ActiveBendDOF)=FEM.modeshape(1:bendModeNo,modeNo);
    dx=min(Xcoord):(max(Xcoord)-min(Xcoord))/100:max(Xcoord);
    dy=min(Ycoord):(max(Ycoord)-min(Ycoord))/100:max(Ycoord);
    [x3,y3]=meshgrid(dx,dy,0);
    z3=griddata(Xcoord,Ycoord,deformUZ,x3,y3,'v4');
    % ZQ=griddata(Xcoord,Ycoord,Zcoord,deformUZ,x3,y3,z3,'linear');
    hf=figure;
    axes1 = axes('Parent',hf,'YTick',[0 0.8 0.801],'XTick',[0 0.6 0.601],...
        'DataAspectRatio',[1 1 1]);
    hold(axes1,'all');
    
    [max_value,max_id]=max(abs(z3(:)));
    set(gcf,'color','w')
    
    scalefactor = 1/(max(abs(z3(:)))/max([max(abs(y3(:))) max(abs(x3(:)))]))/10;
    
    if max(abs(z3(:))) == max((z3(:)))
        sign_max = 1;
    elseif max(abs(z3(:))) == -min((z3(:)))
        sign_max = -1;
    end
    
    
    surf(x3,y3,sign_max*(z3)*scalefactor ,'FaceColor','interp',...
        'LineStyle','none',...
        'EdgeColor','none',...
        'FaceLighting','phong');hold on;
    colorbar('FontSize',16)
    colormap(jet(20));
    %         colorbar('YLim',[-max(unique(z3)),max(unique(z3))]);
    switch Solver
        case 'vibration'
            
            Natural_Freq=frequency(modeNo);
            
            %             title(['Mode shape of Mode ' num2str(modeNo) ', \omega=' num2str(Natural_Freq) 'Hz'],'FontSize',12);
        case 'prestressed_vibr'
            
            Natural_Freq=frequency(modeNo);
            
            %             title(['Prestressed mode shape of Mode ' num2str(modeNo) ', \omega=' num2str(Natural_Freq) 'Hz'],'FontSize',12);
            
        case 'buckling'
            
            %             title(['Buckling Mode ' num2str(modeNo) ', Load Factor \lambda=' num2str(frequency(modeNo))]);%,...
            %             sprintf('\n'),'(\delta=' num2str(delta) ', \gamma=' num2str(gamma) ',\beta=' num2str(beta) ')'],'FontSize',12);
            
            gamma=0;
            beta=0;
            
            %             title(['Mode ' num2str(modeNo) ', Load Factor \lambda=' num2str(frequency(modeNo)),...
            %                 sprintf('\n'),'(ds/hs=' num2str(depthratio(depthNO)) ', \gamma=' num2str(gamma) ',\beta=' num2str(beta) ')'],'FontSize',12);
    end
    for kk=1:size(XXstiffener,1)
        ZZstiffener(kk,:)=ones(1,size(XXstiffener(kk,:),2))*max([Stiffener.height max(abs(z3(:)))*scalefactor]);
        plot3(XXstiffener(kk,:),YYstiffener(kk,:),ZZstiffener(kk,:),'k-','LineWidth',3);hold on;
    end
    hold off;
    xlabel(['Length of the plate, a'],'FontSize',12);
    ylabel(['Width of the plate, b'],'FontSize',12);
    axis image;axis off;
    
    switch Solver
        
        case 'prestressed_vibr'
            
            if lambda_b_ratio>0
                
                export_fig 
                
                filename=['vibrmodeshape_Nxy_depthraio_' num2str(depthratio(depthNO)) '_tensile.png'];
                
                copyfile('export_fig_out.png',filename);
                
                
            elseif lambda_b_ratio<0
                
                export_fig 
                
                filename=['vibrmodeshape_Nxy_depthraio_' num2str(depthratio(depthNO)) '_comp.png'];
                
                copyfile('export_fig_out.png',filename);
                
            end;
            
            
            
        case 'bucklings'
            
            export_fig 
            
            filename=['buckmodeshape_Nxy_depthraio_' num2str(depthratio(depthNO)) '.png'];
                
            copyfile('export_fig_out.png',filename);
    end
    
    
    %     mode_figure_name=['present_stiffenerII_mode_' num2str(plotmodenNo)];
    %     saveas(gcf,[mode_figure_name '.fig'])
    
    
    
    % % % %     %% Contour
    % % % %     %-------------------------------------------------------------
    % % % %     view(2);
    % % % %
    % % % %     %%%%% plot the contour
    % % % %     hg=figure;
    % % % %     axes2 = axes('Parent',hg,'YTick',[-1 1.8 1.801],...
    % % % %         'XTick',[-1 1.6 1.601],...
    % % % %         'DataAspectRatio',[1 1 1],...
    % % % %         'PlotBoxAspectRatio',[1 1.33333333333333 10],...
    % % % %         'LineWidth',2,...
    % % % %         'FontSize',20);
    % % % %     hold(axes2,'all');
    % % % %
    % % % %     set(gcf,'color','w')
    % % % %
    % % % %     %     [max_value,max_id]=max(abs(z3(:)));
    % % % %
    % % % %
    % % % %     contour(x3,y3,z3/z3(max_id),15,'LineWidth',3);hold on;
    % % % %     axis image;colorbar;
    % % % %
    % % % %     for kk=1:size(XXstiffener,1)
    % % % %         ZZstiffener(kk,:)=10*abs(XXstiffener(kk,:));
    % % % %         plot3(XXstiffener(kk,:),YYstiffener(kk,:),ZZstiffener(kk,:),'k-','LineWidth',3);hold on;
    % % % %     end
    % % % %     axis([0 0.601 0 0.801])
    % % % %     hold off;
    % % % %     colormap(jet);axis image;colorbar('FontSize',14);axis image; box on;
    % % % %     Z=zeros(size(FEM.elementNodes,1),size(FEM.nodesCord,2));
    % % % %     set(gcf, 'PaperPosition', [0 0 4 4]);
    % % % %     set(gcf, 'PaperSize', [4 4]);
    % % % %     %     saveas(gcf,['AIAAJ_Previbr_FreeMode_Design' Stiffener.Design '_hsbs' num2str(depthratio(depthNO))],'pdf');
end