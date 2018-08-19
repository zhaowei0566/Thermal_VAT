function FEM = mesh_QUAD8_v2(xmesh,ymesh)

% MESH NAUTRUAL SPACE FOR 8-noded QUAD
% by Wei Zhao, April-2018

xmesh_nodes_1 = 2*xmesh+1; xcords_1 = linspace(-1,1,xmesh_nodes_1);
xmesh_nodes_2 = xmesh+1;   xcords_2 = linspace(-1,1,xmesh_nodes_2);



ymesh_nodes_1 = ymesh+1; ycords_1 = linspace(-1,1,ymesh_nodes_1);
ymesh_nodes_2 = ymesh; ycords_2 = linspace(-1+2/(2*ymesh),1-2/(2*ymesh),ymesh_nodes_2);



nodes_number_1  = xmesh_nodes_1 * ymesh_nodes_1;
nodes_number_2  = xmesh_nodes_2* ymesh;

node_cords_1 = zeros(nodes_number_1,4); %[nodeid, x, y, z]

nodes_number_total = nodes_number_1 + nodes_number_2;

node_cords_all = zeros(nodes_number_total,4);





%%
for ii = 1:ymesh_nodes_1
    
    for jj = 1:xmesh_nodes_1
        
        
  
        label = (ii - 1)* xmesh_nodes_1 + jj + xmesh_nodes_2*(ii-1);
 
        
        node_cords_all(label,:) = [label, xcords_1(jj),ycords_1(ii),0];
        
        
        
        label2 = (ii-1)*xmesh_nodes_1 + jj;
        NODES1_LABEL( label2) = label;
        
        
    end
    
    
end


figure;plot(node_cords_all(:,2),node_cords_all(:,3),'ro');




%%


for ii = 1:ymesh
    
    for jj = 1:xmesh_nodes_2
        
        
        label = (ii - 1)* xmesh_nodes_2 + jj + ii*xmesh_nodes_1;
        
        
        node_cords_all(label,:) = [label, xcords_2(jj),ycords_2(ii),0];
        
        
        label2 = (ii-1)*xmesh_nodes_2 + jj;
        NODES2_LABEL(label2) = label;
        
        
    end
    
    
end

hold on;plot(node_cords_all(:,2),node_cords_all(:,3),'bo');




All_node_cords =  node_cords_all;

All_node_cords(:,2) = (All_node_cords(:,2) + 1)/2;
All_node_cords(:,3) = (All_node_cords(:,3) + 1)/2;
%% FORM ELEMENT



elements_number = xmesh*ymesh;

elements_con = zeros(elements_number,8);


for ii = 1:ymesh
    
    
    for jj = 1:xmesh
        
        
      elem_label  =(ii-1)*xmesh + jj;
      
      
      pnt1_id = (ii-1)*(xmesh*2+1) + 2*(jj-1)+1; label_pnt1  = NODES1_LABEL(pnt1_id);
      
      label_pnt2  = NODES1_LABEL(pnt1_id+2);
      label_pnt5  = NODES1_LABEL(pnt1_id+1);
      
      pnt4_id = (ii)*(xmesh*2+1) + 2*(jj-1)+1; label_pnt4 = NODES1_LABEL(pnt4_id);
      label_pnt3 = NODES1_LABEL(pnt4_id+2);
      label_pnt7 = NODES1_LABEL(pnt4_id+1);
      

      
      
      pnt8_id = (ii-1)*(xmesh+1) + jj;
       label_pnt6 = NODES2_LABEL(pnt8_id+1);
        label_pnt8 = NODES2_LABEL(pnt8_id);
        
    
    elements_con(elem_label,:) = [label_pnt1 label_pnt2 label_pnt3 label_pnt4 ...
                                    label_pnt5 label_pnt6 label_pnt7 label_pnt8];
    end
    
    
   
    
end










FEM.nodesCord = All_node_cords;
FEM.elementNodes = elements_con;

FEM.cbar = '';

FEM.numberElements = elements_number;



