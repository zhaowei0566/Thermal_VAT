function FORCE = applied_force_FEM(FEM,Alpha,N0,Plate)

% % % %% ============= OUT-OF-PLANE FORCE ==============================
% % % FEM.Fz=0;% p is the uniform pressure.
% % % FEM.Mx=0;
% % % FEM.My=0;
% % % %%%
% % % FEM.p= -FEM.Fz;
%         ===================== Transverse Loads ==========================
FORCE=LinearForceMatrix(FEM);
%     pointload=1;


%% ============ IN-PLANE LOADS (N) ==========
% Alpha= 2; % Alpha - load type

RHS_nodes = find(FEM.nodesCord(:,2) == Plate.length);

for iii = 1:length(RHS_nodes)
    
    y = FEM.nodesCord(RHS_nodes(iii),3);
    force_temp(RHS_nodes(iii),1) =  (1-Alpha*y/Plate.width) ;
    
end


% N0= 1e3;%%


if Alpha== 2
    
    
    force_function = @(width_y) N0*(1-Alpha*width_y/Plate.width);
    total_force_applied =  integral(force_function,0,Plate.width/2);
    
    force_temp = -force_temp/sum(force_temp(RHS_nodes(1:(length(RHS_nodes)-1)/2,1)))* total_force_applied; %%25.5
    
    
    
else
    bot_load = 1; top_load = 1-Alpha;
    
    
    force_function = @(width_y) N0*(1-Alpha*width_y/Plate.width);
    total_force_applied =  integral(force_function,0,Plate.width);
    
    
    force_total = sum(linspace(top_load,bot_load,101)); % 101 is the nastran model nodes in right edge
    
    
    force_temp = -force_temp/sum(force_temp) *  total_force_applied ;
    
end


FORCE(FEM.NodeNumber*3+1:FEM.NodeNumber*4)=...
    FORCE(FEM.NodeNumber*3+1:FEM.NodeNumber*4)+ force_temp;