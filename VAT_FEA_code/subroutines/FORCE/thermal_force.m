function TF = thermal_force(Mat, Stru,FEM,Laminate,T0T1,center,physical_length,Temperature)

gDOF=FEM.GDof;

TF = zeros(gDOF,1); % thermal force vector


[GaussWeights,GaussLocations]=gaussQuadrature(FEM.GaussPointBend);


for elem=1:FEM.numberElements
    
    NodeIndices=FEM.elementNodes(elem,:);%% Node NO. for one element
    
    numberNodes=size(FEM.nodeCoordinates,1); %% how many nodes in FEM
    
    % Displacement indicies
    % each node has 5 DOFs
    % w, theta_x, theta_y, u, v
    elementDof=[NodeIndices NodeIndices+numberNodes NodeIndices+2*numberNodes ...
        NodeIndices+3*numberNodes NodeIndices+4*numberNodes];
    
    num_elem_nodes = length(NodeIndices); %% how many nodes for one element;
    
    
    
    N_deltaT = zeros(3,1);
    M_deltaT = zeros(3,1);
    
    % Loop for Gauss Points
    for ee=1:size(GaussWeights,1)
        
        
        GaussPoint=GaussLocations(ee,:);
        xi=GaussPoint(1);eta=GaussPoint(2);
        
        % shape function and derivatives
        [shape,naturalderivatives]=shapefunctionshell(xi,eta,FEM);
        
        % Jacobian and XY chain derivatives  % Jacobian Matrix and Inverse of Jacobian Matrix
        [Jacob,invJacob,XYderivatives]= JacobianFEMshell(naturalderivatives,FEM,NodeIndices);
        
        
        %         centerX = sum(FEM.nodeCoordinates_label( NodeIndices,2))/length( NodeIndices);
        centerX  = shape'*FEM.nodeCoordinates( NodeIndices,1);
        centerY  = shape'*FEM.nodeCoordinates( NodeIndices,2);
        
        
        Amatrix = zeros(3,3);
        Bmatrix = zeros(3,3);
        Dmatrix = zeros(3,3);
        Ashear  = zeros(2,2);
        
        
        
        for layer = 1:size(T0T1,1)
            
            T0 = T0T1(layer,1);
            T1 = T0T1(layer,2);
            
            theta = VAT_fiber_ply_angle_1D(T0,T1,centerY,center(2),physical_length);
            
            
            [AmatrixK,DmatrixK,AshearK,BmatrixK,QthetaK,ThermalExpCoeff]=LaminatedComposite_VAT(Mat,Stru,Laminate,layer,theta);
            
            Amatrix = Amatrix+AmatrixK;
            Bmatrix = Bmatrix+BmatrixK;
            Dmatrix = Dmatrix+DmatrixK;
            Ashear  = Ashear +AshearK;
            
            
        end
        
        Qbar_out_of_plane = QthetaK([1 2 5] , [1 2 5]);
        
        panel_thickness = sum(Laminate.thickness);        
        
        N_deltaT =  N_deltaT + (- Qbar_out_of_plane*Temperature.deltaT_b*  panel_thickness * ThermalExpCoeff);
        
        M_deltaT =  M_deltaT +( - Qbar_out_of_plane*Temperature.k1*panel_thickness^3/12*ThermalExpCoeff);
        
        
        
    end
    
    
    
end