function [ActiveDof,SidesDof]=BCLaminate(FEM)
% Boundary condition for simply supported laminated plate
% J. N. Reddy's Example
type=FEM.BCtype;
%% FIND the SidesDof through Code for Square plate or Rectangle Plate
xx=FEM.nodeCoordinates(:,1);
yy=FEM.nodeCoordinates(:,2);
nodeNum=size(FEM.nodeCoordinates,1);
BCtype=FEM.BCReddyType
%%
switch type
    case 'SSSS' % Simply-supported boundary conditions, w=0;
        switch BCtype
            case 'FSS-1' % Cross-Ply
                % Sym
                % x=0,a; v=w=betaY=0;
                % y=0,b; u=w=betaX=0;
                
                SideNodesNumListY = find(yy==min(FEM.nodeCoordinates(:,2))|...
                    yy==max(FEM.nodeCoordinates(:,2)));
                SideNodesNumListX = find(xx==max(FEM.nodeCoordinates(:,1))|...
                    xx==min(FEM.nodeCoordinates(:,1)));
                
                SidesDofX=[SideNodesNumListX SideNodesNumListX+2*nodeNum SideNodesNumListX+4*nodeNum];
                SidesDofY=[SideNodesNumListY SideNodesNumListY+nodeNum SideNodesNumListY+3*nodeNum ];
                
                SidesDof=[SidesDofX SidesDofY];
                ActiveDof=setdiff(1:FEM.GDof,unique(SidesDof));
                
            case 'FSS-2'
                
                % x=0,a; u=w=betaY=0;
                % y=0,b; v=w=betaX=0;
                SideNodesNumListY = find(yy==min(FEM.nodeCoordinates(:,2))|...
                    yy==max(FEM.nodeCoordinates(:,2)));
                SideNodesNumListX = find(xx==max(FEM.nodeCoordinates(:,1))|...
                    xx==min(FEM.nodeCoordinates(:,1)));
                
                SidesDofX=horzcat(SideNodesNumListX', SideNodesNumListX'+2*nodeNum ,SideNodesNumListX'+3*nodeNum);
                SidesDofY=horzcat(SideNodesNumListY', SideNodesNumListY'+nodeNum, SideNodesNumListY'+4*nodeNum );
                
                SidesDof=horzcat(SidesDofX, SidesDofY);
                ActiveDof=setdiff(1:FEM.GDof,SidesDof);
                
            case 'SSSS-F'
                
                % x=0; u=betaX=0;
                % x=a; v=w=betaY=0;
                % y=0; v=betaY=0;
                % y=b; u=w=betaX=0;
                
                SideNodesNumListY = find(yy==min(FEM.nodeCoordinates(:,2))|...
                    yy==max(FEM.nodeCoordinates(:,2)));
                SideNodesNumListX = find(xx==max(FEM.nodeCoordinates(:,1))|...
                    xx==min(FEM.nodeCoordinates(:,1)));
                
                SidesDofX=[SideNodesNumListX  SideNodesNumListX+3*nodeNum SideNodesNumListX+4*nodeNum];
                SidesDofY=[SideNodesNumListY  SideNodesNumListY+3*nodeNum SideNodesNumListY+4*nodeNum];
                
                SidesDof=[SidesDofX SidesDofY];
                ActiveDof=setdiff(1:FEM.GDof,unique(SidesDof));
                
            case 'SSSS-F3'
                
                % x=0; u=betaX=0;
                % x=a; v=w=betaY=0;
                % y=0; v=betaY=0;
                % y=b; u=w=betaX=0;
                
                SideNodesNumListY = find(yy==min(FEM.nodeCoordinates(:,2))|...
                    yy==max(FEM.nodeCoordinates(:,2)));
                SideNodesNumListX = find(xx==max(FEM.nodeCoordinates(:,1))|...
                    xx==min(FEM.nodeCoordinates(:,1)));
                
                SidesDofX=[SideNodesNumListX  SideNodesNumListX+3*nodeNum SideNodesNumListX+4*nodeNum];
                SidesDofY=[SideNodesNumListY  SideNodesNumListY+3*nodeNum SideNodesNumListY+4*nodeNum];
                
                SidesDof=[SidesDofX SidesDofY];
                ActiveDof=setdiff(1:FEM.GDof/5*3,unique(SidesDof));
                
            case 'SS-2'
                
                % x=0; u=betaX=0;
                % x=a; v=w=betaY=0;
                % y=0; v=betaY=0;
                % y=b; u=w=betaX=0;
                
                SideNodesNumListYmin = find(yy==min(FEM.nodeCoordinates(:,2)));
                SideNodesNumListYmax = find(yy==max(FEM.nodeCoordinates(:,2)));
                
                SideNodesNumListXmin = find(xx==min(FEM.nodeCoordinates(:,1)));
                SideNodesNumListXmax = find(xx==max(FEM.nodeCoordinates(:,1)));
                
                SidesDofXmin=[SideNodesNumListXmin+4*nodeNum SideNodesNumListXmin+nodeNum];
                SidesDofXmax=[SideNodesNumListXmax+3*nodeNum  SideNodesNumListXmax  SideNodesNumListXmax+2*nodeNum ];
                SidesDofYmin=[SideNodesNumListYmin+3*nodeNum SideNodesNumListXmin+2*nodeNum];
                SidesDofYmax=[SideNodesNumListYmax+4*nodeNum  SideNodesNumListYmax  SideNodesNumListYmax+nodeNum ];
                
                SidesDof=unique([SidesDofXmin SidesDofXmax SidesDofYmin SidesDofYmax]);
                
                
                ActiveDof=setdiff(1:FEM.GDof,unique(SidesDof));
                
            case 'SSSS-R'
                
                SideNodesNumListXmin = find(xx==min(FEM.nodeCoordinates(:,1))); % x=0, v=w=betaY=0;
                SidesDofXmin=[SideNodesNumListXmin SideNodesNumListXmin+2*nodeNum SideNodesNumListXmin+4*nodeNum];
                
                SideNodesNumListXmax = find(xx==max(FEM.nodeCoordinates(:,1))); % x=a, v=w=betaY=0;
                SidesDofXmax=[SideNodesNumListXmax SideNodesNumListXmax+2*nodeNum SideNodesNumListXmax+4*nodeNum];
                
                SideNodesNumListYmin = find(yy==min(FEM.nodeCoordinates(:,2))); % y=0, u=w=betaX=0;
                SidesDofYmin=[SideNodesNumListYmin SideNodesNumListYmin+1*nodeNum SideNodesNumListXmin+3*nodeNum];
                
                SideNodesNumListYmax = find(yy==max(FEM.nodeCoordinates(:,2))); % y=b, u=w=betaX=0;
                SidesDofYmax=[SideNodesNumListYmax  SideNodesNumListYmax+1*nodeNum  SideNodesNumListYmax+3*nodeNum];
                
                SidesDof=unique([SidesDofXmin SidesDofXmax SidesDofYmin SidesDofYmax]);
                
                ActiveDof=setdiff(1:FEM.GDof,unique(SidesDof));
                
                
            case 'SSFF-R'
                
                SideNodesNumListXmin = find(xx==min(FEM.nodeCoordinates(:,1))); % x=0, v=w=betaY=0;
                SidesDofXmin=[SideNodesNumListXmin SideNodesNumListXmin+2*nodeNum SideNodesNumListXmin+4*nodeNum];
                
                SideNodesNumListXmax = find(xx==max(FEM.nodeCoordinates(:,1))); % x=a, v=w=betaY=0;
                SidesDofXmax=[SideNodesNumListXmax SideNodesNumListXmax+2*nodeNum SideNodesNumListXmax+4*nodeNum];
                
                SideNodesNumListYmin = find(yy==min(FEM.nodeCoordinates(:,2))); % y=0, u=w=betaX=0;
                SidesDofYmin=[];
                
                SideNodesNumListYmax = find(yy==max(FEM.nodeCoordinates(:,2))); % y=b, u=w=betaX=0;
                SidesDofYmax=[];
                
                SidesDof=unique([SidesDofXmin SidesDofXmax SidesDofYmin SidesDofYmax]);
                
                ActiveDof=setdiff(1:FEM.GDof,unique(SidesDof));
                
            case 'SSFS-R'
                
                SideNodesNumListXmin = find(xx==min(FEM.nodeCoordinates(:,1))); % x=0, v=w=betaY=0;
                SidesDofXmin=[SideNodesNumListXmin SideNodesNumListXmin+2*nodeNum SideNodesNumListXmin+4*nodeNum];
                
                SideNodesNumListXmax = find(xx==max(FEM.nodeCoordinates(:,1))); % x=a, v=w=betaY=0;
                SidesDofXmax=[SideNodesNumListXmax SideNodesNumListXmax+2*nodeNum SideNodesNumListXmax+4*nodeNum];
                
                SideNodesNumListYmin = find(yy==min(FEM.nodeCoordinates(:,2))); % y=0, u=w=betaX=0;
                SidesDofYmin=[];
                
                SideNodesNumListYmax = find(yy==max(FEM.nodeCoordinates(:,2))); % y=b, u=w=betaX=0;
                SidesDofYmax=[SideNodesNumListYmax  SideNodesNumListYmax+1*nodeNum  SideNodesNumListYmax+3*nodeNum];
                
                SidesDof=unique([SidesDofXmin SidesDofXmax SidesDofYmin SidesDofYmax]);
                
                ActiveDof=setdiff(1:FEM.GDof,unique(SidesDof));
                
            case 'SSFC-R'
                
                SideNodesNumListXmin = find(xx==min(FEM.nodeCoordinates(:,1))); % x=0, v=w=betaY=0;
                SidesDofXmin=[SideNodesNumListXmin SideNodesNumListXmin+2*nodeNum SideNodesNumListXmin+4*nodeNum];
                
                SideNodesNumListXmax = find(xx==max(FEM.nodeCoordinates(:,1))); % x=a, v=w=betaY=0;
                SidesDofXmax=[SideNodesNumListXmax SideNodesNumListXmax+2*nodeNum SideNodesNumListXmax+4*nodeNum];
                
                SideNodesNumListYmin = find(yy==min(FEM.nodeCoordinates(:,2))); % y=0, u=w=betaX=0;
                SidesDofYmin=[];
                
                SideNodesNumListYmax = find(yy==max(FEM.nodeCoordinates(:,2))); % y=b, u=w=betaX=0;
                SidesDofYmax=[SideNodesNumListYmax  SideNodesNumListYmax+1*nodeNum ...
                    SideNodesNumListYmax+2*nodeNum SideNodesNumListYmax+3*nodeNum SideNodesNumListYmax+4*nodeNum];
                
                SidesDof=unique([SidesDofXmin SidesDofXmax SidesDofYmin SidesDofYmax]);
                
                ActiveDof=setdiff(1:FEM.GDof,unique(SidesDof));
                
            case 'SSSC-R'
                
                SideNodesNumListXmin = find(xx==min(FEM.nodeCoordinates(:,1))); % x=0, v=w=betaY=0;
                SidesDofXmin=[SideNodesNumListXmin SideNodesNumListXmin+2*nodeNum SideNodesNumListXmin+4*nodeNum];
                
                SideNodesNumListXmax = find(xx==max(FEM.nodeCoordinates(:,1))); % x=a, v=w=betaY=0;
                SidesDofXmax=[SideNodesNumListXmax SideNodesNumListXmax+2*nodeNum SideNodesNumListXmax+4*nodeNum];
                
                SideNodesNumListYmin = find(yy==min(FEM.nodeCoordinates(:,2))); % y=0, u=w=betaX=0;
                SidesDofYmin=[SideNodesNumListYmin SideNodesNumListYmin+1*nodeNum SideNodesNumListXmin+3*nodeNum];
                
                SideNodesNumListYmax = find(yy==max(FEM.nodeCoordinates(:,2))); % y=b, u=w=betaX=0;
                SidesDofYmax=[SideNodesNumListYmax  SideNodesNumListYmax+1*nodeNum ...
                    SideNodesNumListYmax+2*nodeNum SideNodesNumListYmax+3*nodeNum SideNodesNumListYmax+4*nodeNum];
                
                SidesDof=unique([SidesDofXmin SidesDofXmax SidesDofYmin SidesDofYmax]);
                
                ActiveDof=setdiff(1:FEM.GDof,unique(SidesDof));
                
            case 'SSCC-R'
                
                SideNodesNumListXmin = find(xx==min(FEM.nodeCoordinates(:,1))); % x=0, v=w=betaY=0;
                SidesDofXmin=[SideNodesNumListXmin SideNodesNumListXmin+2*nodeNum SideNodesNumListXmin+4*nodeNum];
                
                SideNodesNumListXmax = find(xx==max(FEM.nodeCoordinates(:,1))); % x=a, v=w=betaY=0;
                SidesDofXmax=[SideNodesNumListXmax SideNodesNumListXmax+2*nodeNum SideNodesNumListXmax+4*nodeNum];
                
                SideNodesNumListYmin = find(yy==min(FEM.nodeCoordinates(:,2))); % y=0, u=w=betaX=0;
                SidesDofYmin=[SideNodesNumListYmin SideNodesNumListYmin+1*nodeNum ...
                    SideNodesNumListYmin+2*nodeNum SideNodesNumListYmin+3*nodeNum SideNodesNumListYmin+4*nodeNum];
                
                SideNodesNumListYmax = find(yy==max(FEM.nodeCoordinates(:,2))); % y=b, u=w=betaX=0;
                SidesDofYmax=[SideNodesNumListYmax  SideNodesNumListYmax+1*nodeNum ...
                    SideNodesNumListYmax+2*nodeNum SideNodesNumListYmax+3*nodeNum SideNodesNumListYmax+4*nodeNum];
                
                SidesDof=unique([SidesDofXmin SidesDofXmax SidesDofYmin SidesDofYmax]);
                
                ActiveDof=setdiff(1:FEM.GDof,unique(SidesDof));
        end
        
end