%% FORCE MATRIX
% p -- distributed force
% Mx -- distributed moment about x-axis
% My -- distributed moment about y-axis
function [FORCE]=linearforcebeam(Loads,FEM,Stru,p_type,alpha)

numberNodes=size(FEM.nodeCoordinates,1);

nodesCord=FEM.nodeCoordinates(unique(Loads.Elem(:)),:);

gDOF=numberNodes*2;

FORCE=zeros(gDOF,1);

Nxx=Loads.udl(1);
Nyy=Loads.udl(2);

for e=1:size(nodesCord,1)-1
    
    NodeIndices=Loads.Elem(e,:);
    
    XYZ=FEM.nodeCoordinates(NodeIndices,:);
    
    [GaussWeights,GaussLocations]=gaussQuadrature1D('2');
    % Loop for Gauss Points
    
    for ee=1:size(GaussWeights,1)
        
        GaussPoint=GaussLocations(ee,:);
        
        xi=GaussPoint(1);
        
        % shape function and derivatives
        %         [shape,naturalderivatives]=shapefunctionbeam(xi,FEM);
        Loads.typestiff='CBAR2';
        
        [shape,naturalderivatives,d2Nds2]=shapefunctionbeam(xi,Loads);
        
        detJ=sqrt((naturalderivatives*XYZ(:,1))^2+(naturalderivatives*XYZ(:,2))^2);

        
        %% Nxx perpendicular to y-axis
        
        px = Nxx*InPlane_DistributedLoad(shape'*XYZ(:,2),p_type,alpha,Stru.width);
        
        FORCE(NodeIndices)=FORCE(NodeIndices)+shape*px*detJ*GaussWeights(ee);
        
        
        %% Nyy perpendicular to x-axis
        py = Nyy*InPlane_DistributedLoad(shape'*XYZ(:,1),p_type,alpha,Stru.length);
        
        FORCE(NodeIndices+numberNodes)=FORCE(NodeIndices+numberNodes*1)+...
            +shape*py*detJ*GaussWeights(ee);
        
    end
    
    
    
end
