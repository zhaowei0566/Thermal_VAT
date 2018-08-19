%% FORCE MATRIX
% at edge perpendicuar to x-axis and force is along x-axis
%
%
function Inplane_F=linearforcebeam_X(Loads,FEM,Stru,p_type,alpha)

numberNodes=size(FEM.nodeCoordinates,1);

% nodesCord=FEM.nodeCoordinates(unique(Loads.Elem(:)),:);

Inplane_gDOF=numberNodes*2;

Inplane_F=zeros(Inplane_gDOF,1);

Nxx=Loads.udl(1);
Nxy=Loads.udl(2);


for elem = 1:size(Loads.Elem,1)
    
    NodeIndices=Loads.Elem(elem,:);
    
    XYZ=FEM.nodeCoordinates(NodeIndices,:);
    
    [GaussWeights,GaussLocations]=gaussQuadrature1D('1');
    % Loop for Gauss Points
    
    for gauss_pnt_num=1:size(GaussWeights,1)
        
        GaussPoint=GaussLocations(gauss_pnt_num,:);
        
        xi = GaussPoint(1);
        
        % shape function and derivatives
        %         [shape,naturalderivatives]=shapefunctionbeam(xi,FEM);
        Loads.typestiff='CBAR2';
        
        [shape,naturalderivatives,d2Nds2]=shapefunctionbeam(xi,Loads);
        
        detJ=sqrt((naturalderivatives*XYZ(:,1))^2+(naturalderivatives*XYZ(:,2))^2);

        
        %% Nxx perpendicular 
        
        px = Nxx*InPlane_DistributedLoad(shape'*XYZ(:,2),p_type,alpha,Stru.width);
        
        Inplane_F(NodeIndices)=Inplane_F(NodeIndices)+shape*px*detJ*GaussWeights(gauss_pnt_num);
        
        %% Nxy
        
        py = Nxy*InPlane_DistributedLoad(shape'*XYZ(:,2),p_type,alpha,Stru.length);
        
        Inplane_F(NodeIndices+numberNodes)=Inplane_F(NodeIndices+numberNodes*1)+...
            +shape*py*detJ*GaussWeights(gauss_pnt_num);
        
    end
    

end
