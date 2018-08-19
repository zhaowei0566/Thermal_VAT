function cord_new=Rz(cord,theta,cord0)



x0=cord0(1);
y0 = cord0(2);
z0 = cord0(3);


s= sind(theta);
c = cosd(theta);

T = [c -s 0;
     s c 0;
     0 0 1];
 
 
 vector_old =[cord(1)-x0, cord(2)-y0, cord(3)- z0];
 
 vector_rotated = T*vector_old';
 
 
 cord_new = vector_rotated' + [x0 y0 z0];
 
 
 
 
 



