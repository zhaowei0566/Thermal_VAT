%% FORCE MATRIX
% Fx -- distributed forces  N/m along x-axis
% Fy -- distributed forces N/m along y-axis
function Inplane_F=InplaneForce(FEM,Stru,type,alpha)

xx=FEM.nodeCoordinates(:,1);
yy=FEM.nodeCoordinates(:,2);

% FORCE=zeros(FEM.GDof,1);


UDL=FEM.UDL;

%% findout sides where the applied loads
% SideNodesNumList_LHS = find(xx==min(FEM.nodeCoordinates(:,1)));
SideNodesNumList_RHS = find(xx>=0.99*max(FEM.nodeCoordinates(:,1)));
SideNodesNumList_THS = find(yy>=0.99*max(FEM.nodeCoordinates(:,2)));
% SideNodesNumList_BHS = find(yy==min(FEM.nodeCoordinates(:,2)));


%%
figure(200);hold on;
plot(FEM.nodeCoordinates(SideNodesNumList_RHS,1),FEM.nodeCoordinates(SideNodesNumList_RHS,2),'go');
hold on;
plot(FEM.nodeCoordinates(SideNodesNumList_THS,1),FEM.nodeCoordinates(SideNodesNumList_THS,2),'go');


%%
Nodes_RHS=FEM.nodesCord(SideNodesNumList_RHS,:); % Nxx
Nodes_THS=FEM.nodesCord(SideNodesNumList_THS,:); % Nyy

% applied load nodes coordinates
Loads.Nxx_nodes=Nodes_RHS;
Loads.Nyy_nodes=Nodes_THS;

% element connectivity nodes
Loads.NxxElem=zeros(size(Nodes_RHS,1)-1,2);
Loads.NyyElem=zeros(size(Nodes_THS,1)-1,2);



for ii=1:size(Loads.NyyElem,1)
    
    NxxElem(ii,:)=Nodes_RHS([ii ii+1],1);
    NyyElem(ii,:)=Nodes_THS([ii ii+1],1);
    
end


%% Nxx and Nxy for line perpendicular to x-axis

Loads.udl=UDL(1,:); % line perpendicular to x-axis with force along x-axis N/m


Loads.Elem = NxxElem;

Inplane_F=linearforcebeam_X(Loads,FEM,Stru,type,alpha);

% Inplane_FORCE(FEM.NodeNumber*3+1:FEM.NodeNumber*5)=Inplane_F;

%% Nyx and Nyy for line perpendicular to y-axis

% Loads.type='NYY';


Loads.Elem = NyyElem;

Loads.udl=UDL(2,:); % Nyx and Nyy line perpendicular to y-axis with force along y-axis

% Loads.elementNodes=Loads.NyyElem;

Inplane_F = Inplane_F + linearforcebeam_Y(Loads,FEM,Stru,type,alpha);



