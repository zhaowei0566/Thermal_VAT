function [K,Kelem]=LinearStiffnessLaminatedPlate_VAT_v2_Y(Mat, Stru,FEM,Laminate,T0T1,center,physical_length)
% elementype=1: CTRIA-3
% elementype=2: CTRIA-6
% elementype=3: CQUAD-4
% elementype=4: CQUAD-8
nodedof=FEM.PlateNodeDof;

gDOF=FEM.GDof;

[K]=zeros(gDOF,gDOF);

% Cycle for each element, element stiffness
for elem=1:FEM.numberElements
    NodeIndices=FEM.elementNodes(elem,:);%% Node NO. for one element
    numberNodes=size(FEM.nodeCoordinates,1); %% how many nodes in FEM
    
    % Displacement indicies
    % each node has 5 DOFs
    % w, theta_x, theta_y, u, v
    elementDof=[NodeIndices NodeIndices+numberNodes NodeIndices+2*numberNodes ...
        NodeIndices+3*numberNodes NodeIndices+4*numberNodes];
    num_elem_nodes=length(NodeIndices); %% how many nodes for one element;
    
    Kelem=zeros(num_elem_nodes*5,num_elem_nodes*5);%%% Element Stiffness Matrix;
    
    % ------------Bending Matrix-------------------
    [GaussWeights,GaussLocations]=gaussQuadrature(FEM.GaussPointBend);
    % Loop for Gauss Points
    for ee=1:size(GaussWeights,1)
        
        
        GaussPoint=GaussLocations(ee,:);
        xi=GaussPoint(1);eta=GaussPoint(2);
        
        % shape function and derivatives
        [shape,naturalderivatives]=shapefunctionshell(xi,eta,FEM);
        % Jacobian and XY chain derivatives  % Jacobian Matrix and Inverse of Jacobian Matrix
        [Jacob,invJacob,XYderivatives]=...
            JacobianFEMshell(naturalderivatives,FEM,NodeIndices);
        %%
        
        %         centerX = sum(FEM.nodeCoordinates_label( NodeIndices,2))/length( NodeIndices);
        Element_center_X  = shape'*FEM.nodeCoordinates( NodeIndices,1);
        Element_center_Y  = shape'*FEM.nodeCoordinates( NodeIndices,2);
        
        
        Amatrix = zeros(3,3);
        Bmatrix = zeros(3,3);
        Dmatrix = zeros(3,3);
        Ashear  = zeros(2,2);
        
        
        
        for layer = 1:size(T0T1,1)
            
            T0 = T0T1(layer,1);
            T1 = T0T1(layer,2);
            
            theta = VAT_fiber_ply_angle_1D(T0,T1,Element_center_Y,center,physical_length);
            
            
            [AmatrixK,DmatrixK,AshearK,BmatrixK,QthetaK,ThermalExpCoeff]=LaminatedComposite_VAT(Mat,Stru,Laminate,layer,theta);
            
            Amatrix = Amatrix+AmatrixK;
            Bmatrix = Bmatrix+BmatrixK;
            Dmatrix = Dmatrix+DmatrixK;
            Ashear  = Ashear +AshearK;
            
            
        end
        
%         Bmatrix
%         Amatrix
%         Dmatrix
        %% --------Membrane strain energy-------------
        Bstretch=zeros(3,num_elem_nodes*5);
        Bstretch(1,num_elem_nodes*3+1:num_elem_nodes*4)=XYderivatives(1,:);
        Bstretch(2,num_elem_nodes*4+1:num_elem_nodes*5)=XYderivatives(2,:);
        Bstretch(3,num_elem_nodes*3+1:num_elem_nodes*4)=XYderivatives(2,:);
        Bstretch(3,num_elem_nodes*4+1:num_elem_nodes*5)=XYderivatives(1,:);
        % B MATRIX FOR STRETCH
        
        K(elementDof,elementDof)=K(elementDof,elementDof)+...
            Bstretch'*Amatrix*Bstretch*GaussWeights(ee)*det(Jacob);%% Membrane strain energy
        
        %% -----Coupling extension-bending strain energy-------------
        Bbending=zeros(3,5*num_elem_nodes);
        Bbending(1,num_elem_nodes+1:num_elem_nodes*2)=XYderivatives(1,:);
        Bbending(2,num_elem_nodes*2+1:num_elem_nodes*3)=XYderivatives(2,:);
        Bbending(3,num_elem_nodes+1:num_elem_nodes*2)=XYderivatives(2,:);
        Bbending(3,num_elem_nodes*2+1:num_elem_nodes*3)=XYderivatives(1,:);
        %
        K(elementDof,elementDof)=K(elementDof,elementDof)+...
            Bstretch'*Bmatrix*Bbending*GaussWeights(ee)*det(Jacob);
        K(elementDof,elementDof)=K(elementDof,elementDof)+...
            Bbending'*Bmatrix*Bstretch*GaussWeights(ee)*det(Jacob);
        %% ---------Bending strain energy-----------------------
        K(elementDof,elementDof)=K(elementDof,elementDof)+...
            Bbending'*Dmatrix*Bbending*GaussWeights(ee)*det(Jacob);
    end
    
    %% --------------- For shear matrix----------------------
    [GaussWeights,GaussLocations]=gaussQuadrature(FEM.GaussPointShear);
    
    for ee=1:size(GaussWeights,1)
        
        
        GaussPoint=GaussLocations(ee,:);
        xi=GaussPoint(1);eta=GaussPoint(2);
        
        % shape function and derivatives
        [shape,naturalderivatives]=shapefunctionshell(xi,eta,FEM);
        
        % Jacobian and XY chain derivatives  % Jacobian Matrix and Inverse of Jacobian Matrix
        [Jacob,invJacob,XYderivatives]=...
            JacobianFEMshell(naturalderivatives,FEM,NodeIndices);
        %%
        %         centerX = sum(FEM.nodeCoordinates_label( NodeIndices,2))/length( NodeIndices);
        Element_center_X  = shape'*FEM.nodeCoordinates( NodeIndices,1);
        Element_center_Y  = shape'*FEM.nodeCoordinates( NodeIndices,2);
        
        Amatrix = zeros(3,3);
        Bmatrix = zeros(3,3);
        Dmatrix = zeros(3,3);
        Ashear =zeros(2,2);
        
        for layer = 1:size(T0T1,1)
            
            T0 = T0T1(layer,1);
            T1 = T0T1(layer,2);
            
            theta = VAT_fiber_ply_angle_1D(T0,T1,Element_center_Y,center,physical_length);
            
            
            [AmatrixK,DmatrixK,AshearK,BmatrixK,QthetaK,ThermalExpCoeff]=LaminatedComposite_VAT(Mat,Stru,Laminate,layer,theta);
            
            Amatrix = Amatrix+AmatrixK;
            Bmatrix = Bmatrix+BmatrixK;
            Dmatrix = Dmatrix+DmatrixK;
            Ashear  = Ashear +AshearK;
            
            
        end
        
        % [B] matrix for bending
        Bshear=zeros(2,5*num_elem_nodes);
        Bshear(1,1:num_elem_nodes)=XYderivatives(1,:);
        Bshear(1,1+num_elem_nodes:num_elem_nodes*2)=shape;
        Bshear(2,1:num_elem_nodes)=XYderivatives(2,:);
        Bshear(2,1+num_elem_nodes*2:num_elem_nodes*3)=shape;
        %
        K(elementDof,elementDof)=K(elementDof,elementDof)+ Bshear'*Ashear*Bshear*GaussWeights(ee)*det(Jacob);
        
        Kelem=Kelem+K(elementDof,elementDof);
    end
    
    det(Kelem);
end


