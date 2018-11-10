function FORCE_thermal = calculate_thermal_force(FEM,Laminate,Mat,Stru,T0T1,center,physical_length,Thermal)

gDOF = FEM.GDof;

FORCE_thermal = zeros(gDOF,1);

[GaussWeights,GaussLocations]=gaussQuadrature(FEM.GaussPointBend);

for elem =1:FEM.numberElements
    NodeIndices=FEM.elementNodes(elem,:);
    
    numberNodes=size(FEM.nodeCoordinates,1);
    % Displacement indicies
    elementDof=[NodeIndices NodeIndices+numberNodes NodeIndices+2*numberNodes NodeIndices+3*numberNodes NodeIndices+4*numberNodes];
   
    num_elem_nodes=length(NodeIndices); %% how many nodes for one element;
    
    centerX = sum(FEM.nodeCoordinates_label(NodeIndices,2))/length( NodeIndices);
    
    %% For each laminate layer
    
    for layer=1:size(T0T1,1)
        
        %     [AmatrixK,DmatrixK,AshearK,BmatrixK,QthetaK,ThermalExpCoeff]=...
        %         LaminatedComposite(Mat,Stru,Laminate,layer);
        
        % middle of each layer, from bottom to up
        
        
        Temp_mid = Thermal.Temp_bot;
        k1       = Thermal.k1;
        
        zcoordinate = -Stru.thickness/2 +Laminate.thickness(layer)/2 +(layer-1)*Laminate.thickness(layer);
        
        delta_Temp  =  Temp_mid*( zcoordinate/Stru.thickness*k1 +1);
%         
%         zbot = -Stru.thickness/2  +(layer-1)*Laminate.thickness(layer);
%         ztop = -Stru.thickness/2  + layer*Laminate.thickness(layer);
%         
        
        T0 = T0T1(layer,1);
        T1 = T0T1(layer,2);
        
        theta = VAT_fiber_ply_angle_1D(T0,T1,centerX,center,physical_length);
        
        [AmatrixK,DmatrixK,AshearK,BmatrixK,Qbar_K,ThermalExpCoeff]=LaminatedComposite_VAT(Mat,Stru,Laminate,layer,theta);
        
        
        %% in-plane stress
        thermal_strain=  ThermalExpCoeff*delta_Temp;
        
%         stress0 = Qbar_K([1 2 5],[1 2 5])*strain0;
        
%         thermal_resultant_k = (ztop-zbot)*stress0; % Nxx,Nxy,Nzz, N/m
        
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
            %% --------Membrane strain energy-------------
            Bstretch=zeros(3,num_elem_nodes*5);
            Bstretch(1,num_elem_nodes*3+1:num_elem_nodes*4)=XYderivatives(1,:);
            Bstretch(2,num_elem_nodes*4+1:num_elem_nodes*5)=XYderivatives(2,:);
            Bstretch(3,num_elem_nodes*3+1:num_elem_nodes*4)=XYderivatives(2,:);
            Bstretch(3,num_elem_nodes*4+1:num_elem_nodes*5)=XYderivatives(1,:);
            
            
            FORCE_thermal(elementDof)  =  FORCE_thermal(elementDof) +...
                Bstretch'*AmatrixK*thermal_strain*det(Jacob)*GaussWeights(ee);
            
        end
    end
end