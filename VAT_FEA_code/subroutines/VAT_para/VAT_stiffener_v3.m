function [stiffener_cords,XY_stiffener] = VAT_stiffener_v3(T01,a,b,stiffener_nodes_number)


number_of_nodes = 10001;

x_RHS = linspace(0,a/2,number_of_nodes );

T0 = T01(1)/180*pi;
T1 = T01(2)/180*pi;

if T0 == T1
    
    T1 = T1+1e-2/180*pi;
    
end

% RHS
y_RHS = a/2/(T1-T0).*(-log(cos(T0+2*(T1-T0).*x_RHS/a))+log(cos(T0)));


x_LHS = -fliplr(x_RHS);
y_LHS = -fliplr(y_RHS);

stiffener_cords(1,:) = [x_LHS(1:number_of_nodes-1) x_RHS];
stiffener_cords(2,:) = [y_LHS(1:number_of_nodes-1) y_RHS];



x_label_lower = find(stiffener_cords(1,:)>=-a/2);
x_label_upper = find(stiffener_cords(1,:)<=a/2);

x_labels      = intersect(x_label_lower,x_label_upper);


stiffener_cords = stiffener_cords(:,x_labels);


% y_label_lower = find(stiffener_cords(2,:)>=-b/2);
%
% y_label_upper = find(stiffener_cords(2,:)<=b/2);
%

%% cut from y



y_label_low = find(stiffener_cords(2,:)<-b/2);
y_label_up = find(stiffener_cords(2,:)>b/2);

%%

if ~isempty(y_label_low) && ~isempty(y_label_up)
    
    stiffener_cords = stiffener_cords (:,y_label_up(end):y_label_low(1));
    
end
figure;plot(stiffener_cords(1,:),stiffener_cords(2,:));axis([-a/2 a/2 -b/2 b/2])

%%

X_stiffener = linspace(-max(abs(stiffener_cords(1,:))),max(abs(stiffener_cords(1,:))),stiffener_nodes_number);
Y_stiffener = interp1(stiffener_cords(1,:),stiffener_cords(2,:),X_stiffener, 'spline' );
%% interpolate the node coordinates for the given nodes number

% find which is close to a/2 or b/2

% % % dis_x = max(abs(stiffener_cords(1,:))) - a/2;
% % % dis_y = max(abs(stiffener_cords(2,:))) - b/2;
% % %
% % % if abs(dis_x) > abs(dis_y)
% % %
% % % % %      [min_value,minid] = min(stiffener_cords(1,:));
% % % % %      [max_value,maxid] = max(stiffener_cords(1,:));
% % % % %   if minid >1 && min_value == -a/2
% % % % %
% % % % %       minid = minid - 1;
% % % % %   end
% % % % %   if maxid < size(stiffener_cords,2) && max_value == b/2
% % % % %      maxid = maxid +1;
% % % % %   end
% % % % %
% % % % %
% % %     Y_stiffener = linspace(-max(abs(stiffener_cords(2,:))),max(abs(stiffener_cords(2,:))),stiffener_nodes_number);
% % %
% % %     X_stiffener = interp1(stiffener_cords(2,:),stiffener_cords(1,:),Y_stiffener, 'spline' );
% % %
% % % else
% % %
% % % %     [min_value,minid] = min(stiffener_cords(2,:));
% % % %     [max_value,maxid] = max(stiffener_cords(2,:));
% % % %
% % % %    if minid >1
% % % %       minid = minid - 1;
% % % %
% % % %   end
% % % %   if maxid < size(stiffener_cords,2)
% % % %      maxid = maxid +1;
% % % %   end
% % %
% % %     X_stiffener = linspace(-max(abs(stiffener_cords(1,:))),max(abs(stiffener_cords(1,:))),stiffener_nodes_number);
% % %
% % %     Y_stiffener = interp1(stiffener_cords(1,:),stiffener_cords(2,:),X_stiffener, 'spline' );
% % %
% % % end

XY_stiffener=[X_stiffener;Y_stiffener];

