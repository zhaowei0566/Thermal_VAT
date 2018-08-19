% function MAC = cal_MAC (vec1,vec2)


function mAc=cal_MAC(Phi1,Phi2)



% This function calculates mac between phi1 and phi2
mAc= (abs(Phi1'*Phi2))^2/((Phi1'*Phi1)*(Phi2'*Phi2));

end