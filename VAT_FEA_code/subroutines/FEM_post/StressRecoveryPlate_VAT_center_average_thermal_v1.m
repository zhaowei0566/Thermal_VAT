function [stress_global,strain_global]=StressRecoveryPlate_VAT_center_average_thermal_v1(FEM,...
    Laminate,Mat,Stru,T0T1,center,width,Thermal)
% used to calculate the stress distribution from linear static analysis
% Column 1 - \sigma_{xx}
% Column 2 - \sigma_{yy}
% Column 3 - \tau_{xy}


% u = u0 + z \beta_x
% v = v0 + z \beta_y
% w = w0

% stress recovery is conducted at center of each plate
% z coordinate for each laminate layer is the middle one of the top and
% bottom layer.

disp ( ' ========T H E R M A L      S T R E S S      R E C O V E R Y ==========');
% deformUZ=FEM.displacement(1:FEM.GDof/5);
% deformBx=FEM.displacement(1+FEM.GDof/5:2*FEM.GDof/5);
% deformBy=FEM.displacement(1+2*FEM.GDof/5:3*FEM.GDof/5);
% deformUX=FEM.displacement(1+3*FEM.GDof/5:4*FEM.GDof/5);
% deformUY=FEM.displacement(1+4*FEM.GDof/5:5*FEM.GDof/5);
%


%% For each laminate layer

for layer=1:size(T0T1,1)
    
    %     [AmatrixK,DmatrixK,AshearK,BmatrixK,QthetaK,ThermalExpCoeff]=...
    %         LaminatedComposite(Mat,Stru,Laminate,layer);
    
    % middle of each layer, from bottom to up
    
    
    Temp_mid = Thermal.Temp_bot;
    k1       = Thermal.k1;

    
    
    zcoordinate = -Stru.thickness/2 +Laminate.thickness(layer)/2 +(layer-1)*Laminate.thickness(layer);
    
   delta_Temp  =  Temp_mid*( zcoordinate/Stru.thickness*k1 +1);
    
%     delta_Temp = Temp_bot + k1*(zcoordinate -  (-Stru.thickness/2));
    
    
    for elem=1:FEM.numberElements
        
        NodeIndices=FEM.elementNodes(elem,:);%% Node NO. for one element
        numberNodes=size(FEM.nodeCoordinates,1); %% how many nodes in FEM
        

        
        % Displacement indicies
        % each node has 5 DOFs
        % w, theta_x, theta_y, u, v
        
        
        
        
        
        elementDof=[NodeIndices NodeIndices+numberNodes NodeIndices+2*numberNodes ...
            NodeIndices+3*numberNodes NodeIndices+4*numberNodes];
        
        num_elem_nodes=length(NodeIndices); %% how many nodes for one element;
        
        [GaussWeights,GaussLocations]=gaussQuadrature(FEM.GaussPointBend);
        

        
        NormalStrain=zeros(3,size(GaussWeights,1)); %
        
        %         GaussLocations = [0 0]; GaussWeights = 1;
        
        %%  Loop for Gauss Points for in-plane stress \sigma_xx, \sigma_yy; and \sigma_{xy}
        
        
        for gauss_point_num=1:size(GaussWeights,1)
            
            
            % no need to sum them; strain at these Gaussian Points
            
            GaussPoint=GaussLocations(gauss_point_num,:);
            
            xi= GaussPoint(1);
            eta=GaussPoint(2);
            
            % shape function and derivatives
            [shape,naturalderivatives]=shapefunctionshell(xi,eta,FEM);
            
            % Jacobian and XY chain derivatives  % Jacobian Matrix and Inverse of Jacobian Matrix
            [Jacob,invJacob,XYderivatives] = JacobianFEMshell(naturalderivatives,FEM,NodeIndices);
            
            %% VAT
            
            %         centerX = sum(FEM.nodeCoordinates_label( NodeIndices,2))/length( NodeIndices);
            centerX  = shape'*FEM.nodeCoordinates( NodeIndices,1);
            centerY  = shape'*FEM.nodeCoordinates( NodeIndices,2);
            
            
            
            T0 = T0T1(layer,1);
            T1 = T0T1(layer,2);
            
            theta = VAT_fiber_ply_angle_1D(T0,T1,centerX,center,width);
            
            [AmatrixK,DmatrixK,AshearK,BmatrixK,Qbar_K,ThermalExpCoeff]=LaminatedComposite_VAT(Mat,Stru,Laminate,layer,theta);
            

            
            %% in-plane stress
            strain0(:,gauss_point_num)=  - ThermalExpCoeff * delta_Temp;
            
            stress0(:,gauss_point_num)=Qbar_K([1 2 5],[1 2 5])*strain0(:,gauss_point_num);
            

            
            
        end
        
        
        
        stress_global(elem,:,layer) =  sum(stress0,2)'/gauss_point_num;
        strain_global(elem,:,layer)=   sum(strain0,2)'/gauss_point_num;
        
    end
    
end

