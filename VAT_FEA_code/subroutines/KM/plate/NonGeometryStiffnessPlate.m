function [KG]=NonGeometryStiffnessPlate(FEM,Laminate,Stru)
% Gometric stiffness obtained from calculated stress [sigma_xx,sigma_yy,tau_xy]
gDOF=FEM.GDof;
KG=zeros(gDOF,gDOF);




switch FEM.typeplate
    
    case 'CQUAD8'
        
        
        %-----Geometry Stiffness due to Tranverse deformation------
        % Cycle for each element, element stiffness,
        
        for e=1:FEM.numberElements
            
            for layer=1:size(Laminate.thickness,2)
                
                
                Zcoordinate1=-Stru.thickness/2+(layer-1)*Laminate.thickness(layer);
                Zcoordinate2=-Stru.thickness/2+layer*Laminate.thickness(layer);
                
                stress0=[FEM.stress(e,1,layer),FEM.stress(e,3,layer);
                         FEM.stress(e,3,layer),FEM.stress(e,2,layer)];
                
                
                NodeIndices=FEM.elementNodes(e,:);%% Node NO. for one element
                numberNodes=size(FEM.nodeCoordinates,1);
                % Displacement indicies
                % each node has 3 DOFs
                elementDof=[NodeIndices NodeIndices+numberNodes NodeIndices+2*numberNodes...
                    NodeIndices+3*numberNodes NodeIndices+4*numberNodes];
                num_elem_nodes=length(NodeIndices); %% how many nodes for one element;
                
                [GaussWeights,GaussLocations]=gaussQuadrature('2by2');
                % Loop for Gauss Points
                for ee=1:size(GaussWeights,1)
                    
                    GaussPoint=GaussLocations(ee,:);
                    xi=GaussPoint(1);eta=GaussPoint(2);
                    
                    % shape function and derivatives
                    [shape,naturalderivatives]=shapefunctionshell(xi,eta,FEM);
                    
                    % Jacobian and XY chain derivatives  % Jacobian Matrix and Inverse of Jacobian Matrix
                    [Jacob,invJacob,XYderivatives]=...
                        JacobianFEMshell(naturalderivatives,FEM,NodeIndices);
                    
                    Lsbending=zeros(2,8);
                    Lsbending(1,1:8)=XYderivatives(1,:);
                    Lsbending(2,1:8)=XYderivatives(2,:);
                    
                    KG(NodeIndices,NodeIndices)=KG(NodeIndices,NodeIndices)+...
                        (Zcoordinate2-Zcoordinate1)*Lsbending'*stress0*Lsbending*GaussWeights(ee)*det(Jacob);
                end
                % ---------Geometry STIFFNESS due to inplane residual stress--
                
                [GaussWeights,GaussLocations]=gaussQuadrature('2by2');
                
                % Loop for Gauss Points
                for ee=1:size(GaussWeights,1)
                    
                    GaussPoint=GaussLocations(ee,:);
                    xi=GaussPoint(1);eta=GaussPoint(2);
                    
                    % shape function and derivatives
                    [shape,naturalderivatives]=shapefunctionshell(xi,eta,FEM);
                    
                    % Jacobian and XY chain derivatives  % Jacobian Matrix and Inverse of Jacobian Matrix
                    [Jacob,invJacob,XYderivatives]=...
                        JacobianFEMshell(naturalderivatives,FEM,NodeIndices);
                    
                    Lshear1=zeros(2,8);
                    Lshear1(1,1:8)=XYderivatives(1,:);Lshear1(2,1:8)=XYderivatives(2,:);
                    KG(NodeIndices+numberNodes,NodeIndices+numberNodes)=KG(NodeIndices+numberNodes,NodeIndices+numberNodes)+...
                        (Zcoordinate2^3-Zcoordinate1^3)/3*Lshear1'*stress0*Lshear1*GaussWeights(ee)*det(Jacob);
                    %
                    
                    Lshear2=zeros(2,8);
                    Lshear2(1,1:8)=XYderivatives(1,:);Lshear2(2,1:8)=XYderivatives(2,:);
                    KG(NodeIndices+numberNodes*2,NodeIndices+numberNodes*2)=KG(NodeIndices+numberNodes*2,NodeIndices+numberNodes*2)+...
                    (Zcoordinate2^3-Zcoordinate1^3)/3*Lshear2'*stress0*Lshear2*GaussWeights(ee)*det(Jacob);
                 %      Laminate.thickness(layer)^3/12*Lshear2'*stress0*Lshear2*GaussWeights(ee)*det(Jacob);
                end
            end
            
        end
end












