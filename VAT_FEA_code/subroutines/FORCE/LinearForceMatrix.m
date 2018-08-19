%% FORCE MATRIX
% p -- distributed transverse force
% Mx -- distributed moment about x-axis
% My -- distributed moment about y-axis
function [FORCE]=LinearForceMatrix(FEM)
gDOF=FEM.GDof;

[FORCE]=zeros(gDOF,1);

% Fx=FEM.Fx;
% Fy=FEM.Fy;
p=FEM.p;
Mx=FEM.Mx;
My=FEM.My;

[GaussWeights,GaussLocations]=gaussQuadrature(FEM.GaussPointBend);

for e=1:FEM.numberElements
    NodeIndices=FEM.elementNodes(e,:);
    
    numberNodes=size(FEM.nodeCoordinates,1);
    % Displacement indicies
    %     elementDof=[NodeIndices NodeIndices+numberNodes NodeIndices+2*numberNodes NodeIndices+3*numberNodes NodeIndices+4*numberNodes];
    num_elem_nodes=length(NodeIndices); %% how many nodes for one element;
    
    % Loop for Gauss Points
    for ee=1:size(GaussWeights,1)
        
        GaussPoint=GaussLocations(ee,:);
        xi=GaussPoint(1);eta=GaussPoint(2);
        
        % shape function and derivatives
        [shape,naturalderivatives]=shapefunctionshell(xi,eta,FEM);
        
        % Jacobian and XY chain derivatives  % Jacobian Matrix and Inverse of Jacobian Matrix
        [Jacob,invJacob,XYderivatives]=...
            JacobianFEMshell(naturalderivatives,FEM,NodeIndices);
        %
        %% Distributed Forces & Moments; Fz Mx My Fx Fy
        FORCE(NodeIndices)=FORCE(NodeIndices)+shape*p*det(Jacob)*GaussWeights(ee);
        
        FORCE(NodeIndices+numberNodes)=FORCE(NodeIndices+numberNodes)+...
            shape*Mx*det(Jacob)*GaussWeights(ee);
        
        FORCE(NodeIndices+2*numberNodes)=FORCE(NodeIndices+2*numberNodes)+...
            shape*My*det(Jacob)*GaussWeights(ee);
        
% % %         FORCE(NodeIndices+3*numberNodes)=FORCE(NodeIndices+3*numberNodes)+...
% % %             shape*Fx*det(Jacob)*GaussWeights(ee);
% % %         
% % %         FORCE(NodeIndices+4*numberNodes)=FORCE(NodeIndices+4*numberNodes)+...
% % %             shape*Fy*det(Jacob)*GaussWeights(ee);
        
        %% FORCE - Moment_x
        % FORCE(NodeIndices+numberNodes)=
        %% FORCE - Moment_y
        % FORCE(NodeIndices+2*numberNodes)=
    end
end
