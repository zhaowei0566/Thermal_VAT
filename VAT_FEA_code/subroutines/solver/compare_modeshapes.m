function [value,id,mAc]= compare_modeshapes(target_modes,Ktotal,KGplate,ActiveDof)




[Vp,Dp]=eigs(Ktotal(ActiveDof,ActiveDof),-KGplate(ActiveDof,ActiveDof),10,'sm');

[DDP,modeNop]=sort((diag(Dp)));

VVsortP=Vp(:,modeNop);

mineiglabel=find(min(abs(DDP))==-DDP);
if isempty(mineiglabel)==1
    mineiglabel=find(min(abs(DDP))==DDP);
end
plotmodenNo=mineiglabel;

buckling_mode_shape = VVsortP(:,plotmodenNo);



for ii = 1:size(target_modes,2)

    target_mode  = target_modes(:,ii);
    
    mAc(ii) = cal_MAC(buckling_mode_shape,target_mode);
    
end


[value,id] = max(mAc);






