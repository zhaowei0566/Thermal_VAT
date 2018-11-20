function [Mass,Melem]=LinearMassLaminatedPlate(FEM,Mat,Stru)

% Mass matrix for isotropic and laminated composite panel

gDOF=FEM.GDof;

rho=Mat.density;

Mass=zeros(gDOF,gDOF);

% Cycle for each element, element stiffness
for e=1:FEM.numberElements
    
    NodeIndices=FEM.elementNodes(e,:);%% Node NO. for one element
    numberNodes=size(FEM.nodeCoordinates,1); %% how many nodes in FEM
    
    % Displacement indicies
    % each node has 5 DOFs
    % w, theta_x, theta_y, u, v
    
    elementDof=[NodeIndices NodeIndices+numberNodes NodeIndices+2*numberNodes ...
        NodeIndices+3*numberNodes NodeIndices+4*numberNodes];
    
    num_elem_nodes=length(NodeIndices); %% how many nodes for one element;
    
    Melem=zeros(num_elem_nodes*5,num_elem_nodes*5);%%% Element Stiffness Matrix;
    
    % ------------Bending Matrix-------------------
    [GaussWeights,GaussLocations]=gaussQuadrature(FEM.GaussPointBend);
    % Loop for Gauss Points
    for ee=1:size(GaussWeights,1)
        
        GaussPoint=GaussLocations(ee,:);
        
        xi=GaussPoint(1);eta=GaussPoint(2);
        
        % shape function and derivatives
        
        [shape,naturalderivatives]=shapefunctionshell(xi,eta,FEM);
        
        % Jacobian and XY chain derivatives
        % Jacobian Matrix and Inverse of Jacobian Matrix
        [Jacob,invJacob,XYderivatives]=...
            JacobianFEMshell(naturalderivatives,FEM,NodeIndices);
        
        %% --------Structural Mass Matrix-------------
        
        NmassShape=zeros(5,num_elem_nodes*5);
        
        NmassShape(1,1:num_elem_nodes)=shape';
        
        NmassShape(2,1+num_elem_nodes:num_elem_nodes*2)=shape';
        
        NmassShape(3,1+num_elem_nodes*2:num_elem_nodes*3)=shape';
        
        NmassShape(4,1+num_elem_nodes*3:num_elem_nodes*4)=shape';
        
        NmassShape(5,1+num_elem_nodes*4:num_elem_nodes*5)=shape';
        
        
        
        t=Stru.thickness;
        
%         mp=rho*[t, 0, 0, 0, 0;
%             0, t, 0, 0, 0;
%             0, 0, t, 0, 0;
%             0, 0, 0, t^3/12, 0;
%             0, 0, 0, 0, t^3/12];
        
        mp=rho*[t,0,0,0,0;
                0,t^3/12,0,0,0;
                0,0,t^3/12,0,0;
                0,0,0,t,0;
                0,0,0,0,t];
        
        
        Mass(elementDof,elementDof)=Mass(elementDof,elementDof)+...
            NmassShape'*mp*NmassShape*GaussWeights(ee)*det(Jacob);
        
        Melem=Melem+NmassShape'*mp*NmassShape*GaussWeights(ee)*det(Jacob);
    end    
    
end


