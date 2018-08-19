function FORCE=DiscreteForceMatrix(FEM,Stru,pointload)

%% Y-Stiffener

FORCE=zeros(FEM.GDof,1);

XX=FEM.nodesCord(:,2);
YY=FEM.nodesCord(:,3);


for nodeID=1:FEM.NodeNumber
    
    
    %     if XX(nodeID)<Stru.length/2
    %         FORCE(nodeID)=pointload;
    %
    %     elseif XX(nodeID)>Stru.length/2
    %
    %
    %         FORCE(nodeID)=-pointload;
    %     end
    
    if XX(nodeID)==Stru.length/2 && YY(nodeID)==Stru.width/2
        
        %         FORCE(nodeID+FEM.GDof/5)=pointload;
        FORCE(nodeID+FEM.GDof/5*1)=pointload;
    end
    
end