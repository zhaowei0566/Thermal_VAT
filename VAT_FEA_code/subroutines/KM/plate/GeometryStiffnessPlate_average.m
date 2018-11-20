function KG=GeometryStiffnessPlate_average(FEM,Stru,Laminate)
gDOF=FEM.GDof;
KG=zeros(gDOF,gDOF);




switch FEM.typeplate
        
    case 'CQUAD8'

% % %         stress0=[FEM.stress0(1,1),FEM.stress0(1,2);
% % %             FEM.stress0(2,1),FEM.stress0(2,2)];
        
        %-----Geometry Stiffness due to Tranverse deformation------
        % Cycle for each element, element stiffness,
        for elem=1:FEM.numberElements
            
            % ================
            
            stress_x = (Laminate.thickness * squeeze(FEM.stress(elem,1,:)))/sum(Laminate.thickness);
            
            stress_y = (Laminate.thickness * squeeze(FEM.stress(elem,2,:)))/sum(Laminate.thickness);
            stress_xy = (Laminate.thickness * squeeze(FEM.stress(elem,3,:)))/sum(Laminate.thickness);
            
            stress0 = [stress_x, stress_xy;
                           stress_xy, stress_y];
            
            % ===============
            
            NodeIndices=FEM.elementNodes(elem,:);%% Node NO. for one element
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
                Lsbending(1,1:8)=XYderivatives(1,:);Lsbending(2,1:8)=XYderivatives(2,:);
                KG(NodeIndices,NodeIndices)=KG(NodeIndices,NodeIndices)+...
                    Stru.thickness*Lsbending'*stress0*Lsbending*GaussWeights(ee)*det(Jacob);
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
                    Stru.thickness^3/12*Lshear1'*stress0*Lshear1*GaussWeights(ee)*det(Jacob);
                %
                
                Lshear2=zeros(2,8);
                Lshear2(1,1:8)=XYderivatives(1,:);Lshear2(2,1:8)=XYderivatives(2,:);
                KG(NodeIndices+numberNodes*2,NodeIndices+numberNodes*2)=KG(NodeIndices+numberNodes*2,NodeIndices+numberNodes*2)+...
                    Stru.thickness^3/12*Lshear2'*stress0*Lshear2*GaussWeights(ee)*det(Jacob);
            end
            
            
            
        end
 
end













