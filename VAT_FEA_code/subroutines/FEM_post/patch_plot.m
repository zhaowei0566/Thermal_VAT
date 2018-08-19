function patch_plot(t,p,FID,componment)
% plot triangle or quadlateral mesh based on points and element
% connectivity
% Dec-17.
% input:
%   t: element connectivity - just element id; no label, no pid
%   p: nodes for those elements - include label and coordinates
%   FID: figure ID for plot in one figure.
% The input t or p probably has very arbitry label, it is better to sort
% them in correct format for patch


% 4-node element or 3-node element

% remove all zeros in t



%%
switch componment
    
    case 'skin'
        
        facecolor = [0.800000011920929 0.800000011920929 0.800000011920929];
        edgecolor = [1 0 0];
    case 'spar'
        facecolor = [1 1 1];
        edgecolor = [0 1 0];
    case 'rib'
        
         facecolor = [1 1 0];
         edgecolor = [0 0 0];
         
    case 'mesh'
        
        facecolor = [0.313725501298904 0.313725501298904 0.313725501298904];
        edgecolor = [ 1 0 0];
end

%%
jj=1;

for t_num_temp = 1:size(t,1)
    
    if all(t(t_num_temp ,:))
        
        t_updated(jj,:) = t(t_num_temp,:);
        
        jj=jj+1;
        
    end
end

% t = t_updated;

%% sort nodal label.

[sorted_nodal_label,sorted_id] = sort(p(:,1)); % p should have unique label here.

p_updated_only_cordinates = p(sorted_id,2:4);


%%
element_nodes_num = size(t,2);

elements_num = size(t,1);




if element_nodes_num == 4 || element_nodes_num == 8
    
    % needs to know the mixed triangle elements
    
    ctria3_t_index = find(t(:,4) == 0);
    
    if isempty(ctria3_t_index) == 0
        
        t_ctria3 = t( ctria3_t_index,1:3);
        
        for elem_num_temp = 1:size(t_ctria3,1)
            
            for node_num_temp = 1:3
                
                % find nodal label in sorted_nodal_label index
                
                index_temp = find(t_ctria3(elem_num_temp,node_num_temp) == sorted_nodal_label);
                
                t_ctria3(elem_num_temp,node_num_temp) = index_temp;%sorted_id(index_temp);
            end
            
        end
        
        figure(FID);hold on;
        patch('Faces',t_ctria3,'Vertices',p_updated_only_cordinates,'FaceColor',facecolor,'FaceAlpha',.3,'EdgeColor',edgecolor);
        
        
    end
    
    %% cquad4
    
    left_index = setdiff(1:elements_num,ctria3_t_index);
    t_cquad4 = t(left_index ,1:4);
    
    
    for elem_num_temp = 1:size(t_cquad4,1)
        
        for node_num_temp = 1:4
            
            % find nodal label in sorted_nodal_label index
            
            index_temp = find(t_cquad4(elem_num_temp,node_num_temp) == sorted_nodal_label);
            
            t_cquad4(elem_num_temp,node_num_temp) = index_temp;% sorted_id(index_temp);
        end
        
    end
    
    
    figure(FID);hold on;
    patch('Faces',t_cquad4,'Vertices',p_updated_only_cordinates,'FaceColor',facecolor,'FaceAlpha',.3,'EdgeColor',edgecolor);
    
    
    %% all triangle elements
elseif element_nodes_num == 3
    
    t_temp = t;
    
    for elem_num_temp = 1:size(t_temp,1)
        
        for node_num_temp = 1:3
            
            % find nodal label in sorted_nodal_label index
            
            index_temp = find(t_temp(elem_num_temp,node_num_temp) == sorted_nodal_label);
            
            t_temp(elem_num_temp,node_num_temp) = index_temp;%sorted_id(index_temp);
            
        end
        
        
    end
   
    %% plot
    figure(FID);hold on;
    patch('Faces',t_temp,'Vertices',p_updated_only_cordinates,'FaceColor',facecolor,'FaceAlpha',.3,'EdgeColor',edgecolor);
    
end


