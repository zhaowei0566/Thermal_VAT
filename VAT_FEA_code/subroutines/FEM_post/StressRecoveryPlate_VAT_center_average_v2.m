function [stress_global,strain_global]=StressRecoveryPlate_VAT_center_average_v2(FEM,Laminate,Mat,Stru,T0T1,center,width)
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

disp ( ' ======== S T R E S S      R E C O V E R Y ==========');
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
    
    % middle of each layer
    zcoordinate = -Stru.thickness/2 +Laminate.thickness(layer)/2 +(layer-1)*Laminate.thickness(layer);
    
    
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
        Curvature0=zeros(3,size(GaussWeights,1));
        ShearStrain=zeros(2,size(GaussWeights,1));
        
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
            
  
            
            %% --------------------- In-plane strain at middle surface------------------
            % w, beta_x, beta_y, u, v
            % epsilon_x, epsilon_y, gamma_{xy}
            Bstretch=zeros(3,num_elem_nodes*5);
            
            Bstretch(1,num_elem_nodes*3+1:num_elem_nodes*4)=XYderivatives(1,:);
            Bstretch(2,num_elem_nodes*4+1:num_elem_nodes*5)=XYderivatives(2,:);
            
            Bstretch(3,num_elem_nodes*3+1:num_elem_nodes*4)=XYderivatives(2,:);
            Bstretch(3,num_elem_nodes*4+1:num_elem_nodes*5)=XYderivatives(1,:);
            
            %             NormalStrain=NormalStrain+Bstretch*FEM.displacement(elementDof)';
            NormalStrain(:,gauss_point_num)=Bstretch*FEM.displacement(elementDof)'; % in-plane strain due to in-plane displacement
            %% -----Coupling extension-bending strain energy at middle surface----------
            % kappa_x, kappa_y and kappa_{xy}
            Bbending=zeros(3,5*num_elem_nodes);
            Bbending(1,num_elem_nodes+1:num_elem_nodes*2)   = XYderivatives(1,:);
            Bbending(2,num_elem_nodes*2+1:num_elem_nodes*3) = XYderivatives(2,:);
            
            Bbending(3,num_elem_nodes+1:num_elem_nodes*2)  =  XYderivatives(2,:);
            Bbending(3,num_elem_nodes*2+1:num_elem_nodes*3) = XYderivatives(1,:);
            
            Curvature0(:,gauss_point_num) = Bbending*FEM.displacement(elementDof)';
            
            
            
            %% in-plane stress
            strain0(:,gauss_point_num)=[NormalStrain(1,gauss_point_num)+zcoordinate*Curvature0(1,gauss_point_num),...
                                        NormalStrain(2,gauss_point_num)+zcoordinate*Curvature0(2,gauss_point_num), ...
                                        NormalStrain(3,gauss_point_num)+zcoordinate*Curvature0(3,gauss_point_num)]';
            
            stress0(:,gauss_point_num)=Qbar_K([1 2 5],[1 2 5])*strain0(:,gauss_point_num);
            
            

            
            
            %% Shear strain is not considered since it is out-of-plane strain
            
            
        end
        
        
%                     if elem == 1
%                 
%                 strain0
%                 
%                 stress0
%                 
%                 
%                 
%             end
        %         size(stress0)
        
        stress_global(elem,:,layer) =  sum(stress0,2)'/gauss_point_num;
        strain_global(elem,:,layer)=   sum(strain0,2)'/gauss_point_num;
        
    end
    
end

