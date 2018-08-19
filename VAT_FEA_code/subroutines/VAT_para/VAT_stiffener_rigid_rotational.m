function [stiffener_cords,XY_stiffener] = VAT_stiffener_rigid_rotational(T01,a,b,stiffener_nodes_number,RigidTheta)


number_of_nodes = 10001;

x_RHS = linspace(0,a/2,number_of_nodes );

T0 = T01(1)/180*pi;
T1 = T01(2)/180*pi;

if T0 == T1
    
    T1 = T1+1e-2/180*pi;
    
end

% RHS w.r.t x'
y_RHS = a/2/(T1-T0).*(-log(cos(T0+2*(T1-T0).*x_RHS/a))+log(cos(T0)));


x_LHS = -fliplr(x_RHS);
y_LHS = -fliplr(y_RHS);

stiffener_cords(1,:) = [x_LHS(1:number_of_nodes-1) x_RHS];
stiffener_cords(2,:) = [y_LHS(1:number_of_nodes-1) y_RHS];


figure;plot(stiffener_cords(1,:),stiffener_cords(2,:))
%%


stiffener_cords_rotated = stiffener_cords;

for ii = 1:size(stiffener_cords,2)
    
    stiffener_cords_temp = stiffener_cords(:,ii);
    stiffener_cords_rotated_temp =Rz([stiffener_cords_temp' 0] ,RigidTheta,[0 0 0]);
    
    stiffener_cords_rotated(:,ii) = stiffener_cords_rotated_temp(1:2)';
    
end

hold on;plot(stiffener_cords_rotated(1,:),stiffener_cords_rotated(2,:));

%%


x_label_lower = find(stiffener_cords(1,:)>=-a/2);
x_label_upper = find(stiffener_cords(1,:)<=a/2);
x_labels      = intersect(x_label_lower,x_label_upper);



y_label_lower = find(stiffener_cords(2,:)>=-b/2);
y_label_upper = find(stiffener_cords(2,:)<=b/2);
y_labels      = intersect(y_label_lower,y_label_upper);


final_labels  = intersect(x_labels,y_labels);


stiffener_cords = stiffener_cords(:,final_labels);


figure;plot(stiffener_cords(1,:),stiffener_cords(2,:))
%% interpolate the node coordinates for the given nodes number

% find which is close to a/2 or b/2

dis_x = max(abs(stiffener_cords(1,:))) - a/2;
dis_y = max(abs(stiffener_cords(2,:))) - b/2;

if abs(dis_x) > abs(dis_y)
    
    [value,minid] = min(stiffener_cords(1,:));
    [value,maxid] = max(stiffener_cords(1,:));
    if minid >1
        minid = minid - 1;
    end
    if maxid < size(stiffener_cords,2)
        maxid = maxid +1;
    end
    Y_stiffener = linspace(-max(abs(stiffener_cords(2,minid:maxid))),max(abs(stiffener_cords(2,minid:maxid))),stiffener_nodes_number);
    
    X_stiffener = interp1(stiffener_cords(2,minid:maxid),stiffener_cords(1,minid:maxid),Y_stiffener, 'spline' );
    
else
    
    [value,minid] = min(stiffener_cords(2,:));
    [value,maxid] = max(stiffener_cords(2,:));
    
    if minid >1
        minid = minid - 1;
        
    end
    if maxid < size(stiffener_cords,2)
        maxid = maxid +1;
    end
    
    X_stiffener = linspace(-max(abs(stiffener_cords(1,maxid:minid))),max(abs(stiffener_cords(1,maxid:minid))),stiffener_nodes_number);
    
    Y_stiffener = interp1(stiffener_cords(1,maxid:minid),stiffener_cords(2,maxid:minid),X_stiffener, 'spline' );
    
end

XY_stiffener=[X_stiffener;Y_stiffener];

