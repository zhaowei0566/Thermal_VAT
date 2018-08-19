% postprocess for solution from stiffenedPlate
modeshapez=zeros(1,FEM.GDof);
modeshapez(1:FEM.GDof)=FEM.displacement;

deformUZ=FEM.displacement(1:FEM.GDof/5);
deformBx=FEM.displacement(1+FEM.GDof/5:2*FEM.GDof/5);
deformBy=FEM.displacement(1+2*FEM.GDof/5:3*FEM.GDof/5);



deformUX=FEM.displacement(1+3*FEM.GDof/5:4*FEM.GDof/5);
deformUY=FEM.displacement(1+4*FEM.GDof/5:5*FEM.GDof/5);

% max(deformUZ)
% max(deformUX)
% max(deformUY)
% max(deformBx)
% max(deformBy)



% s1 = find(FEM.nodeCoordinates_label(:,3) == Stru.width/2)
% 
% temp2 = FEM.nodeCoordinates_label(s1,:)
% 
% s2 = find(temp2(:,2) == Stru.length/2);
% 
% temp3 = temp2(s2,:)
% 
% deformUX(temp3(1))
% deformUY(temp3(1))





% Spline deflection and angle
dx=min(Xcoord):(max(Xcoord)-min(Xcoord))/100:max(Xcoord);
dy=min(Ycoord):(max(Ycoord)-min(Ycoord))/100:max(Ycoord);
[x3,y3]=meshgrid(dx,dy,0);

x1 =  griddata(Xcoord,Ycoord,deformUX,x3,y3,'v4');
y2 =  griddata(Xcoord,Ycoord,deformUY,x3,y3,'v4');
z3 =  griddata(Xcoord,Ycoord,deformUZ,x3,y3,'v4');
BetaX=griddata(Xcoord,Ycoord,deformBx,x3,y3,'v4');
BetaY=griddata(Xcoord,Ycoord,deformBy,x3,y3,'v4');

figure
surf(x3,y3,x1,'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong');colormap(jet);view(2);axis image;
title('In-plane displacement, UX');colorbar

figure
surf(x3,y3,y2,'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong');colormap(jet);view(2);axis image;
title('In-plane displacement, UY');colorbar

figure
surf(x3,y3,z3,'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong');colormap(jet);view(2);axis image;
title('In-plane displacement, UZ');colorbar

figure
surf(x3,y3,BetaX,'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong');colormap(jet);view(2);axis image;
title('In-plane displacement, Beta_X');colorbar


figure
surf(x3,y3,BetaY,'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong');colormap(jet);view(2);axis image;
title('In-plane displacement, Beta_Y');colorbar

% griddata(Xcoord,Ycoord,deformUX,Stru.length/2,Stru.width/2,'v4')

% griddata(Xcoord,Ycoord,deformUY,Stru.length/2,Stru.width/2,'v4')

% %% Plot in-plane deflections
% figure(111);hold on;
% % plot(Xcoord,Ycoord,'k.');hold on;
% title('Nodes of Stiffened Plate','FontSize',15);axis image;hold on;
% SF=0.1*Stru.length/max(abs([deformUX deformUY deformUZ]));
% 
% 
% FEM.nodeCoordinates_label_deformed = FEM.nodeCoordinates_label;
% FEM.nodeCoordinates_label_deformed(:,2) = Xcoord'+SF*deformUX;
% FEM.nodeCoordinates_label_deformed(:,3) = Ycoord'+SF*deformUY;
% 
% % deformed body
% patch_plot(FEM.elementNodes,FEM.nodeCoordinates_label_deformed,200,'mesh')
%%  Transverse Deflection Plot
% % % figure
% % % surf(x3,y3,z3,'FaceColor','interp',...
% % %     'EdgeColor','none',...
% % %     'FaceLighting','phong');
% % % if ebar==0
% % %     title(['Transverse Deflection of stiffened plate under uniform load, Concentric'],'FontSize',14);
% % % else
% % %     title(['Transverse Deflection of stiffened plate under uniform load, Eccentric'],'FontSize',14);
% % % end
% % % 
% % % xlabel('Length of the plate, a /m');
% % % ylabel('Width of the plate, b/m');
% % % zlabel('Transverse deflection /m');
% % % colorbar('location','SouthOutside');
% % % colormap('jet')
% % % view(2);
% % % axis image;

disp('The Largest Displacement UX is:');
[value,id] = max(abs(deformUX));
deformUX(id)

disp('The Largest Displacement UY is:');
[value,id] = max(abs(deformUY));
deformUY(id)



disp('The Largest Transverse Displacement is:');
[value,id] = max(abs(deformUZ));
deformUZ(id)



disp('The Largest Transverse Displacement is:');
[value,id] = max(abs(deformUZ));
deformUZ(id)


disp('The Largest Rotations BetaX is:');
[value,id] = max(abs(deformBx));
deformBx(id)

disp('The Largest Rotations BetaY is:');
[value,id] = max(abs(deformBy));
deformBy(id)
%% Plot Transverse Deflection along center line
% % XXcord=x3(51,:); YYcord=y3(:,51); ZZcordx=z3(:,51); %%x=0.06;
% % ZZcordy=z3(51,:);
% % figure(201);hold on
% % plot(XXcord, ZZcordx,'k-','LineWidth',2);
% % % legend('stiffener, present analysis')
% % xlabel(['Distance along centerline(m)'],'FontSize',12);
% % ylabel(['Deflection(m)'],'FontSize',12);
% % hold off;axis normal;
% % 
% % 
% % figure(100);hold on;
% % plot(YYcord,ZZcordy,'r-','LineWidth',2);
% % % legend('Present,x=a/2', 'FontSize',12);
% % xlabel(['Distance along centerline(m)'],'FontSize',12);
% % ylabel(['Deflection(m)'],'FontSize',12);
% % hold off;axis normal;




