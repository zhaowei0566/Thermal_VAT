
function [x_rotated,y_rotated] = Rz_group(cord,Phi,cord0)

x =cord(:,1);
y = cord(:,2);

for ii = 1:length(x)

    cord_new=Rz([x(ii) y(ii) 0],Phi/pi*180,cord0);
    
    x_rotated(ii) = cord_new(1);
    y_rotated(ii) = cord_new(2);

end