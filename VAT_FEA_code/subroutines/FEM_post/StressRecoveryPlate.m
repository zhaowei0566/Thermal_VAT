function [stress,strain]=StressRecoveryPlate(FEM,Laminate,Mat,Stru)
% used to calculate the stress distribution from linear static analysis
% Column 1 - \sigma_{xx}
% Column 2 - \sigma_{yy}
% Column 3 - \tau_{xy}


% u = u0 + z \beta_x
% v = v0 + z \beta_y
% w = w0



deformUZ=FEM.displacement(1:FEM.GDof/5);
deformBx=FEM.displacement(1+FEM.GDof/5:2*FEM.GDof/5);
deformBy=FEM.displacement(1+2*FEM.GDof/5:3*FEM.GDof/5);
deformUX=FEM.displacement(1+3*FEM.GDof/5:4*FEM.GDof/5);
deformUY=FEM.displacement(1+4*FEM.GDof/5:5*FEM.GDof/5);



%% For each laminate layer

for layer=1:length(Laminate.thickness)
    
    [AmatrixK,DmatrixK,AshearK,BmatrixK,Qbar_T,ThermalExpCoeff]=...
        LaminatedComposite(Mat,Stru,Laminate,layer);
    
    zcoordinate=-Stru.thickness/2+(layer-1)*Laminate.thickness(layer);
    
    

    for elem=1:FEM.numberElements
        NodeIndices=FEM.elementNodes(elem,:);%% Node NO. for one element
        numberNodes=size(FEM.nodeCoordinates,1); %% how many nodes in FEM
        
        NormalStrain=zeros(3,1); % 
        Curvature0=zeros(3,1);
        ShearStrain=zeros(2,1);
        
        % Displacement indicies
        % each node has 5 DOFs
        % w, theta_x, theta_y, u, v
        elementDof=[NodeIndices NodeIndices+numberNodes NodeIndices+2*numberNodes ...
            NodeIndices+3*numberNodes NodeIndices+4*numberNodes];
        
        num_elem_nodes=length(NodeIndices); %% how many nodes for one element;

        [GaussWeights,GaussLocations]=gaussQuadrature(FEM.GaussPointBend);
        
        %%  Loop for Gauss Points for in-plane stress \sigma_xx, \sigma_yy; and \sigma_{xy}
        for ee=1:size(GaussWeights,1)
            
            % no need to sum them; strain at these Gaussian Points
            
            GaussPoint=GaussLocations(ee,:);
            
            xi=GaussPoint(1);
            eta=GaussPoint(2);
            
            % shape function and derivatives
            [shape,naturalderivatives]=shapefunctionshell(xi,eta,FEM);
            
            % Jacobian and XY chain derivatives  % Jacobian Matrix and Inverse of Jacobian Matrix
            [Jacob,invJacob,XYderivatives] = JacobianFEMshell(naturalderivatives,FEM,NodeIndices);
            %% --------------------- In-plane strain at middle surface------------------
            % w, beta_x, beta_y, u, v
            % epsilon_x, epsilon_y, gamma_{xy}
            Bstretch=zeros(3,num_elem_nodes*5);
            Bstretch(1,num_elem_nodes*3+1:num_elem_nodes*4)=XYderivatives(1,:);
            Bstretch(2,num_elem_nodes*4+1:num_elem_nodes*5)=XYderivatives(2,:);
            Bstretch(3,num_elem_nodes*3+1:num_elem_nodes*4)=XYderivatives(2,:);
            Bstretch(3,num_elem_nodes*4+1:num_elem_nodes*5)=XYderivatives(1,:);
            
%             NormalStrain=NormalStrain+Bstretch*FEM.displacement(elementDof)';
            NormalStrain(:,ee)=Bstretch*FEM.displacement(elementDof)'; % in-plane strain due to in-plane displacement
            %% -----Coupling extension-bending strain energy at middle surface----------
            % kappa_x, kappa_y and kappa_{xy}
            Bbending=zeros(3,5*num_elem_nodes);
            Bbending(1,num_elem_nodes+1:num_elem_nodes*2)  =XYderivatives(1,:);
            Bbending(2,num_elem_nodes*2+1:num_elem_nodes*3)=XYderivatives(2,:);
            Bbending(3,num_elem_nodes+1:num_elem_nodes*2)  =XYderivatives(2,:);
            Bbending(3,num_elem_nodes*2+1:num_elem_nodes*3)=XYderivatives(1,:);
            
            Curvature0(:,ee) = Bbending*FEM.displacement(elementDof)';


            %% Shear strain is not considered since it is out-of-plane strain
            
            
        end
        
        NormalStrain = sum(NormalStrain,2)/ee;
         Curvature0  = sum(Curvature0,2)/ee ;
        
        strain0=[NormalStrain(1)+zcoordinate*Curvature0(1),...
                 NormalStrain(2)+zcoordinate*Curvature0(2), ...
                 NormalStrain(3)+zcoordinate*Curvature0(3)]';
        
        stress(elem,1:3,layer)=Qbar_T([1 2 5],[1 2 5])*strain0;
        
        strain(elem,1:3,layer)=strain0';
    end
    
end

