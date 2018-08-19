function [All_node_cords,  elements_con] = mesh_QUAD8(xmesh,ymesh)

% MESH NAUTRUAL SPACE FOR 8-noded QUAD


xmesh_nodes_1 = 2*xmesh+1; xcords_1 = linspace(-1,1,xmesh_nodes_1);
xmesh_nodes_2 = xmesh+1;   xcords_2 = linspace(-1,1,xmesh_nodes_2);



ymesh_nodes_1 = ymesh+1; ycords_1 = linspace(-1,1,ymesh_nodes_1);
ymesh_nodes_2 = ymesh; ycords_2 = linspace(-1+2/(2*ymesh),1-2/(2*ymesh),ymesh_nodes_2);



nodes_number_1 = xmesh_nodes_1 * ymesh_nodes_1;
nodes_number_2  = xmesh_nodes_2* ymesh;

node_cords_1 = zeros(nodes_number_1,4); %[nodeid, x, y, z]



%%
for ii = 1:ymesh_nodes_1
    
    for jj = 1:xmesh_nodes_1
        
        label = (ii - 1)* xmesh_nodes_1 + jj;
        
        
        node_cords_1(label,:) = [label, xcords_1(jj),ycords_1(ii),0];
        
        
    end
    
    
end


figure;plot(node_cords_1(:,2),node_cords_1(:,3),'ro');

%%

node_cords_2 = zeros(nodes_number_2,4);

for ii = 1:ymesh
    
    for jj = 1:xmesh_nodes_2
        
        label = (ii - 1)* xmesh_nodes_2 + jj + nodes_number_1;
        
        
        node_cords_2(label,:) = [label, xcords_2(jj),ycords_2(ii),0];
        
        
    end
    
    
end

hold on;plot(node_cords_2(:,2),node_cords_2(:,3),'bo');




All_node_cords =  node_cords_2;

All_node_cords(1:nodes_number_1,:) = node_cords_1;


All_node_cords(:,2) = All_node_cords(:,2) + 1;
All_node_cords(:,3) = All_node_cords(:,3) + 1;
%% FORM ELEMENT



elements_number = xmesh*ymesh;

elements_con = zeros(elements_number,8);


for ii = 1:ymesh
    
    
    for jj = 1:xmesh
        
        
      elem_label  =(ii-1)*xmesh + jj;
        
      element_label = (ii-1)*(2*xmesh+1) + 2*(jj-1)+1;
      
      
      element_label_above = (ii)*(2*xmesh+1) + 2*(jj-1)+1;
      
      
      mid_label_id = (ii-1)*(xmesh+1) + jj +  nodes_number_1;
      
      
      elements_con( elem_label,:) =  [ element_label  element_label+2  element_label_above+2  element_label_above ,...
                                        element_label+1    mid_label_id+1  element_label_above+1 mid_label_id ];
        
        
    
    
    end
    
    
   
    
end














