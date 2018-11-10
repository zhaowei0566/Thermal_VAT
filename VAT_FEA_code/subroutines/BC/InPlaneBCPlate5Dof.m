%% Take off the column in MATRIX with known displacements
% Or you can offer the nodes No. to get the SidesDof
function [ActiveDof,Constrained_Dof]=InPlaneBCPlate5Dof(FEM)

type=FEM.BCtype;
%% FIND the SidesDof through Code for Square plate or Rectangle Plate
xx=FEM.nodeCoordinates(:,1);
yy=FEM.nodeCoordinates(:,2);
nodeNum=size(FEM.nodeCoordinates,1);

%% Offer ActiveDOF manually
% w, betax, betay, u, v
%%

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

switch type
    case 'SSSS-NXX-Thermal'
        
        % Unaxial stress - Nxx
        % left side, u=v=w=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % right side,u=v=w=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_RHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_RHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        
        % top side, u=v=w=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_THS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % bottom side, u=v=w=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_BHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof= unique([SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS]);
        ActiveDof=setdiff(1:FEM.GDof,Constrained_Dof);
        
        
        
    case 'SSSS-NXX'
        
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
        SideNodesNum = SideNodesNum_BHS;
        SideNodesNumList = SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum) = SideNodesNumList;%% w is constrainted
        SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2) = SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof= unique([SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS]);
        ActiveDof=setdiff(1:FEM.GDof,Constrained_Dof);
        
        %         FEM.GDof/5*3+
    case 'SSSS-NXX1'
        
        % Unaxial stress - Nxx
        % left side, u=v=w=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % right side, w=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_RHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % top side, w=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % bottom side, w=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof= unique([SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS]);
        ActiveDof=setdiff(1:FEM.GDof,Constrained_Dof);
        
        
    case 'SSSS-disp-Coburn'
        
        % Unaxial stress - Nxx
        % left side, u=v=w=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % right side, w=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        
        % top side, w=v=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % bottom side, w=v=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof= unique([SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS]);
        ActiveDof=setdiff(1:FEM.GDof,Constrained_Dof);
        
        
    case 'SSSS-disp-Coburn-paper'
        
        % Unaxial stress - Nxx
        % left side, u=v=w=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        %         SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        %         SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % right side, w=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        %         SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        
        % top side, w=v=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        %         SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % bottom side, w=v=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        %         SideNodesNumList=SideNodesNumList_BHS;
        %         SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        %
        SideNodesNumList_MID = find(xx==0.5*max(FEM.nodeCoordinates(:,1)))';
        
        SidesDof_MID =  SideNodesNumList_MID;
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof= unique([SidesDof_LHS  SidesDof_THS SidesDof_BHS  1:FEM.GDof/5*3  ]);
        ActiveDof=setdiff(FEM.GDof/5*3+1:FEM.GDof,Constrained_Dof);
        
        
        
    case 'SSSS-disp-Coburn-v2'
        
        % Unaxial stress - Nxx
        % left side, u=v=w=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        %         SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % right side, w=0
        %         SideNodesNum_RHS=length(SideNodesNumList_RHS);
        %         SideNodesNum=SideNodesNum_RHS;
        %         SideNodesNumList=SideNodesNumList_RHS;
        %         SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        
        % top side, w=v=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        %         SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % bottom side, w=v=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        %         SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof= unique([SidesDof_LHS  SidesDof_THS SidesDof_BHS]);
        ActiveDof=setdiff(1:FEM.GDof,Constrained_Dof);
        
        
    case 'SSSS-NYY'
        % Unaxial stress - NYY
        % left side, u=w=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        %SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % right side, u=w=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_RHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% v is constrainted
        
        % top side, u=w=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% v is constrainted
        
        % bottom side, u=v=w=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_BHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof=unique([SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS]);
        ActiveDof=setdiff(1:FEM.GDof,Constrained_Dof);
        
    case 'SSSS-NYY1'
        % Unaxial stress - NYY
        % left side, w=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        %SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % right side, w=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_RHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% v is constrainted
        
        % top side, w=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% v is constrainted
        
        % bottom side, u=v=w=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_BHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof=unique([SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS]);
        ActiveDof=setdiff(1:FEM.GDof,Constrained_Dof);
        
    case 'SSSS-NXXNYY'
        
        % axial stress - NYY and UXX
        % left side, u=w=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        %SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % right side, w=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_RHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% v is constrainted
        
        % top side, w=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% v is constrainted
        
        % bottom side, v=w=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof=unique([SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS]);
        ActiveDof=setdiff(1:FEM.GDof,Constrained_Dof);
        
        
    case 'SSSS-NXY'
        
        cornernodes=intersect(SideNodesNumList_LHS,SideNodesNumList_BHS);
        % axial stress - NYY and UXX
        % left side, w=0 V=0
        
        
        % left side
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% OPTIONAL -u is constrainted
        %         SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; % V=0
        
        
        
        % right side, w=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_RHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum;
        
        
        % top side, w=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        
        
        % bottom side, w=0, U=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% u is constrainted
        %         SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; % V=0 OPTIMAL
        
        % corner nodes u=v=w=0
        SideNodesNum_corner=length(cornernodes);
        SideNodesNum=SideNodesNum_corner;
        SideNodesNumList=cornernodes;
        SidesDof_Corner(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_Corner(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_Corner(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof=unique([SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS SidesDof_Corner ]);
        ActiveDof=setdiff(1:FEM.GDof,Constrained_Dof);
        
        
    case 'SSSS-NXY2' % not like NASTRAN case
        
        cornernodes=intersect(SideNodesNumList_LHS,SideNodesNumList_BHS);
        % axial stress - NYY and UXX
        % left side, w=0 V=0
        
        
        % left side
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum = SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+4*nodeNum; %% OPTIONAL -u is constrainted
        %         SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; % V=0
        
        
        
        % right side, w=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_RHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% OPTIONAL -u is constrainted
        
        
        % top side, w=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_THS(1+SideNodesNum:2*SideNodesNum)=SideNodesNumList + 4*nodeNum;%% v is constrainted
        
        
        % bottom side, w=0, U=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        %                 SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+3*nodeNum; % V=0 OPTIMAL
        
        % corner nodes u=v=w=0
        SideNodesNum_corner=length(cornernodes);
        SideNodesNum=SideNodesNum_corner;
        SideNodesNumList=cornernodes;
        SidesDof_Corner(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_Corner(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_Corner(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof=unique([SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS SidesDof_Corner ]);
        ActiveDof=setdiff(1:FEM.GDof,Constrained_Dof);
        
    case 'Zhu-Gu'
        % Unaxial stress - NYY
        % left side, x=0; w=betax=v=0
        SideNodesNum_LHS=length(SideNodesNumList_LHS);
        SideNodesNum=SideNodesNum_LHS;
        SideNodesNumList=SideNodesNumList_LHS;
        SidesDof_LHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %  SidesDof_LHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        %SidesDof_LHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % right side, w=0
        SideNodesNum_RHS=length(SideNodesNumList_RHS);
        SideNodesNum=SideNodesNum_RHS;
        SideNodesNumList=SideNodesNumList_RHS;
        SidesDof_RHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_RHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% v is constrainted
        
        % top side, w=0
        SideNodesNum_THS=length(SideNodesNumList_THS);
        SideNodesNum=SideNodesNum_THS;
        SideNodesNumList=SideNodesNumList_THS;
        SidesDof_THS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        %         SidesDof_THS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% v is constrainted
        
        % bottom side, u=v=w=0
        SideNodesNum_BHS=length(SideNodesNumList_BHS);
        SideNodesNum=SideNodesNum_BHS;
        SideNodesNumList=SideNodesNumList_BHS;
        SidesDof_BHS(1:SideNodesNum)=SideNodesNumList;%% w is constrainted
        SidesDof_BHS(1*SideNodesNum+1:SideNodesNum*2)=SideNodesNumList+3*nodeNum; %% u is constrainted
        SidesDof_BHS(2*SideNodesNum+1:SideNodesNum*3)=SideNodesNumList+4*nodeNum; %% v is constrainted
        
        % ActiveDof=setdiff(1:FEM.GDof/5*3,SidesDof);
        Constrained_Dof= unique([SidesDof_LHS SidesDof_RHS SidesDof_THS SidesDof_BHS]);
        ActiveDof=setdiff(1:FEM.GDof,Constrained_Dof);
        
end

