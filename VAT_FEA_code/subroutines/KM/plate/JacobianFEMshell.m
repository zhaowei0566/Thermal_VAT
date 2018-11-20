function [Jacob,invJacob,XYderivatives]=...
    JacobianFEMshell(naturalderivatives,FEM,NodeIndices)


elementype=FEM.typeplate;


switch elementype
    
    case 'CQUAD4'
        Coordinates=FEM.nodeCoordinates;
        
        Jacob=naturalderivatives*Coordinates(NodeIndices,:);
        invJacob=inv(Jacob);
        XYderivatives=invJacob*naturalderivatives;
        
        
    case 'CQUAD8'
        Coordinates=FEM.nodeCoordinates;
        
        Jacob=naturalderivatives*Coordinates(NodeIndices,:);
        invJacob=inv(Jacob);
        XYderivatives=invJacob*naturalderivatives;
       
   
end
