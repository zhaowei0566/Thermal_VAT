%% Take off the column in MATRIX with known displacements
% Or you can offer the nodes No. to get the SidesDof
function [ActiveDof,SidesDof]=EssentialBCPlate5Dof(FEM,Stru)
type=FEM.BCtype;
%% FIND the SidesDof through Code for Square plate or Rectangle Plate
xx=FEM.nodeCoordinates(:,1);
yy=FEM.nodeCoordinates(:,2);
nodeNum=size(FEM.nodeCoordinates,1);



SideNodesNumList_LHS = find(xx==min(FEM.nodeCoordinates(:,1)))';
figure(200);hold on; plot(xx(SideNodesNumList_LHS ),yy(SideNodesNumList_LHS ),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
SideNodesNumList_RHS = find(xx>=0.99*max(FEM.nodeCoordinates(:,1)))';
figure(200);hold on; plot(xx(SideNodesNumList_RHS ),yy(SideNodesNumList_RHS ),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
SideNodesNumList_THS = find(yy>=0.99*max(FEM.nodeCoordinates(:,2)))';
figure(200);hold on; plot(xx(SideNodesNumList_THS ),yy(SideNodesNumList_THS ),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
SideNodesNumList_BHS = find(yy==min(FEM.nodeCoordinates(:,2)))';
figure(200);hold on; plot(xx(SideNodesNumList_BHS ),yy(SideNodesNumList_BHS ),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);



%% Offer ActiveDOF manually
% w, betax, betay, u, v
%%
switch type
    case 'SSSS' % Simply-supported boundary conditions, w=0;
        SideNodesNumList = find(...
            xx==max(FEM.nodeCoordinates(:,1))|...
            xx==min(FEM.nodeCoordinates(:,1))|...
            yy==min(FEM.nodeCoordinates(:,2))|...
            yy==max(FEM.nodeCoordinates(:,2)));
        SideNodesNum=length(SideNodesNumList);
        SidesDof(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        %         ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        ActiveDof=setdiff(1:FEM.GDof,SidesDof);
        
        
    case 'SSSS-inplane'
        
        SideNodesNumList_LHS = find(xx==min(FEM.nodeCoordinates(:,1)))';
        SideNodesNumList_RHS = find(xx>=0.99*max(FEM.nodeCoordinates(:,1)))';
        
        SideNodesNumList_THS = find(yy>=0.99*max(FEM.nodeCoordinates(:,2)))';
        SideNodesNumList_BHS = find(yy==min(FEM.nodeCoordinates(:,2)))';
        
        
        
        % Unaxial stress - Nxx
        % left side, u=v=w=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % right side, v=w=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_RHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % top side, v=w=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % bottom side, v=w=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        SidesDof=[SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS];
        ActiveDof=setdiff(1:FEM.GDof,SidesDof);
        
        %         FEM.GDof/5*3+
        
        
%     case 'SSSS-3' % Simply-supported boundary conditions, w=0;
%         SideNodesNumList = find(...
%             xx>=max(FEM.nodeCoordinates(:,1))-eps|...
%             xx==min(FEM.nodeCoordinates(:,1))|...
%             yy==min(FEM.nodeCoordinates(:,2))|...
%             yy>=max(FEM.nodeCoordinates(:,2))-eps);
%         SideNodesNum=length(SideNodesNumList);
%         
%         
%         
%         SidesDof(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
%         %                 SidesDof(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
%         %                 SidesDof(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
%         ActiveDof=setdiff(1:FEM.GDof,SidesDof);
%         %         ActiveDof=setdiff(1:FEM.GDof,SidesDof);
%         
        
    case 'CCCC'
        SideNodesNumList = find(xx==max(FEM.nodeCoordinates(:,1))|...
            xx==min(FEM.nodeCoordinates(:,1))|...
            yy==min(FEM.nodeCoordinates(:,2))|...
            yy==max(FEM.nodeCoordinates(:,2)));
        SideNodesNum=length(SideNodesNumList);
        
        SidesDof(1:SideNodesNum)=SideNodesNumList;
        SidesDof(SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+nodeNum;
        SidesDof(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+2*nodeNum;
        SidesDof(3*SideNodesNum+1:SideNodesNum*4)=SideNodesNumList+3*nodeNum;
        SidesDof(4*SideNodesNum+1:SideNodesNum*5)=SideNodesNumList+4*nodeNum;
        
        ActiveDof=setdiff(1:FEM.GDof,SidesDof);
        
    case 'CFFF-3'
        SideNodesNumList = find(yy == min(FEM.nodeCoordinates(:,2)));
        SideNodesNum=length(SideNodesNumList);
        SidesDof(1:SideNodesNum)=SideNodesNumList;
        SidesDof(SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+nodeNum;
        SidesDof(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+2*nodeNum;
        SidesDof(3*SideNodesNum+1:SideNodesNum*4)=SideNodesNumList+3*nodeNum;
        SidesDof(4*SideNodesNum+1:SideNodesNum*5)=SideNodesNumList+4*nodeNum;
        ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        
        
    case 'CSCS' %% Modification needed
        SideNodesNumListC = find(xx==max(FEM.nodeCoordinates(:,1))|...
            xx==min(FEM.nodeCoordinates(:,1)));
        SideNodesNumListS=find(yy==min(FEM.nodeCoordinates(:,2))|...
            yy==max(FEM.nodeCoordinates(:,2)));
        SideNodesNumC=length(SideNodesNumListC);
        SideNodesNumS=length(SideNodesNumListS);
        
        SidesDof(1:SideNodesNumC)=SideNodesNumListC;
        SidesDof(SideNodesNumC+1:SideNodesNumC*2)=SideNodesNumListC+nodeNum;
        SidesDof(2*SideNodesNumC+1:SideNodesNumC*3)=SideNodesNumListC+2*nodeNum;
        SidesDof(3*SideNodesNumC+1:SideNodesNumC*4)=SideNodesNumListC+3*nodeNum;
        SidesDof(4*SideNodesNumC+1:SideNodesNumC*5)=SideNodesNumListC+4*nodeNum;
        
        SidesDof(SideNodesNumC*5+1:SideNodesNumC*5+SideNodesNumS)=SideNodesNumListS;
        
        SidesDof(SideNodesNumC*5+SideNodesNumS+1:SideNodesNumC*5+SideNodesNumS*2)=SideNodesNumListS+3*nodeNum;
        SidesDof(SideNodesNumC*5+SideNodesNumS*2+1:SideNodesNumC*5+SideNodesNumS*3)=SideNodesNumListS+4*nodeNum;
        
        ActiveDof=setdiff(1:FEM.GDof,SidesDof);
        
        
    case 'Zhu-Gu' % Simply-supported boundary conditions, w=0;
%         SideNodesNumList_LHS = find(xx==min(FEM.nodeCoordinates(:,1)))';
%         SideNodesNumList_RHS = find(xx>=0.99*max(FEM.nodeCoordinates(:,1)))';
%         
%         SideNodesNumList_THS = find(yy>=0.99*max(FEM.nodeCoordinates(:,2)))';
%         SideNodesNumList_BHS = find(yy==min(FEM.nodeCoordinates(:,2)))';
        
        % left side, w=phi_x=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+2*nodeNum; %% phi_x=0
        
        % right side, w=phi_x=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_RHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+2*nodeNum; %% phi_x=0
        
        % top side, w=phi_y=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+1*nodeNum; %% phi_y=0
        
        % bottom side, w=phyi_y=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+1*nodeNum; %% phi_y=0
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        SidesDof=[SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS];
        ActiveDof=setdiff(1:FEM.GDof/3,SidesDof);
        
        
    case 'SSSS-disp-Coburn'
        
        % Unaxial stress - Nxx
        % left side, u=v=w=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+ nodeNum; %% BETAy is constrainted
        %         SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % right side, w=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_RHS(1+SideNodesNum:2*SideNodesNum)=SideNodesNumList + nodeNum;%% betaY is constrainted
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof= unique([SidesDof_LHS SidesDof_RHS]);
        ActiveDof=setdiff(1:FEM.GDof/5*3,Constrained_Dof);
        
        SidesDof = Constrained_Dof;
        
        case 'SSSS-3'
            
      % left side, w=phi_x=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
%         SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+2*nodeNum; %% phi_x=0
        
        % right side, w=phi_x=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
%         SidesDof_RHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+2*nodeNum; %% phi_x=0
        
        % top side, w=phi_y=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
%         SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+1*nodeNum; %% phi_y=0
        
        % bottom side, w=phyi_y=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
%         SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+1*nodeNum; %% phi_y=0
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        SidesDof=unique([SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS]);
        
        ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        
        
                case 'SSSS-5'
            
      % left side, w=phi_x=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
%         SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+2*nodeNum; %% phi_x=0
        
        % right side, w=phi_x=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
%         SidesDof_RHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+2*nodeNum; %% phi_x=0
        
        % top side, w=phi_y=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
%         SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+1*nodeNum; %% phi_y=0
        
        % bottom side, w=phyi_y=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
%         SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+1*nodeNum; %% phi_y=0
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        SidesDof=unique([SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS]);
        
        ActiveDof=setdiff(1:FEM.GDof,SidesDof);
        
%         SidesDof = Constrained_Dof;
            
    case 'SSSS-disp-Coburn-5dof'
        
        % Unaxial stress - Nxx
        % left side, u=v=w=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+ nodeNum; %% BETAy is constrainted
        SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+3*nodeNum; %% u is constrainted.
        SidesDof_LHS(3*SideNodesNum+1:SideNodesNum*4)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % right side, w=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_RHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+ nodeNum; %% BETAy is constrainted
        SidesDof_RHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+3*nodeNum; %% u is constrainted.
        SidesDof_RHS(3*SideNodesNum+1:SideNodesNum*4)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % top side, w=v=0
        %         SideNodesNum_THS=length(SideNodesNumList_THS);
        %         SideNodesNum=SideNodesNum_THS;
        %         SideNodesNumList=SideNodesNumList_THS;
        %         SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        %
        % bottom side, w=v=0
        %         SideNodesNum_BHS=length(SideNodesNumList_BHS);
        %         SideNodesNum=SideNodesNum_BHS;
        %         SideNodesNumList=SideNodesNumList_BHS;
        %         SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof= unique([SidesDof_LHS SidesDof_RHS]);
        ActiveDof=setdiff(1:FEM.GDof,Constrained_Dof);
        
        SidesDof = Constrained_Dof;
    case 'others'
        
        disp('Input by yourself');
        
        
        
end
