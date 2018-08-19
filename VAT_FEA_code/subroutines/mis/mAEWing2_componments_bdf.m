function mAEWing2_componments_dat(groups,design_folder,design_file_name)


groups_bdf_number=size(groups,2);


for groups_num=1:groups_bdf_number
    
    file_folder_path=[design_folder '\' design_file_name '_' groups{groups_num} '.bdf'];
    
    bdf_node_elem=read_bdf(file_folder_path);
    
    
    fid102=fopen([design_folder '\' design_file_name '_' groups{groups_num} '.dat'],'wt');
    
    for lineNumPCL=1:size(bdf_node_elem,2)
        
        
        fprintf(fid102,'%s\n',[bdf_node_elem{lineNumPCL}]);
        
    end
    
    fclose(fid102);

    
    
end

fclose('all');