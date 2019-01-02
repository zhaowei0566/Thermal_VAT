% The code is mainly used for conducting thermal buckling analysis of VAT
% laminated plate by using example in Duran's paper, 
% Change the material name in Line 104 for different material cases.
%
% Contact: Wei Zhao (weizhao@vt.edu), May, 2018

clear all;warning off;format long;
close all;

% ==== add subroutines ======
%
addpath([pwd filesep 'subroutines' filesep 'KM']);
addpath([pwd filesep 'subroutines' filesep 'KM' filesep 'plate']);
%
addpath([pwd filesep 'subroutines' filesep 'BC']);
addpath([pwd filesep 'subroutines' filesep 'FORCE']);
addpath([pwd filesep 'subroutines' filesep 'mis']);
addpath([pwd filesep 'subroutines' filesep 'FEM_post']);

addpath([pwd filesep 'subroutines' filesep 'FEM_post' filesep 'export_fig']);
%
%
addpath([pwd filesep 'subroutines' filesep 'response']);
addpath([pwd filesep 'subroutines' filesep 'VAT_para']);
%
addpath([pwd filesep 'subroutines' filesep 'mesh']);


global Mat;
global FEM;
global Stru;
global Stiffener;
global Laminate;

global Thermal
global Plate


%% Geometry
Plate.width = 0.15;

Plate.length =0.15;

%% === Meshing plate =====

FEM = mesh_QUAD8_v2(24,24);


%% =============== PANEL GEOMETRY================================
%%%==============================================================
FEM.nodesCord(:,2)=FEM.nodesCord(:,2)/max(FEM.nodesCord(:,2))*Plate.length;
FEM.nodesCord(:,3)=FEM.nodesCord(:,3)/max(FEM.nodesCord(:,3))*Plate.width;
%---------Coordinates of node------------
Xcoord=FEM.nodesCord(:,2);
Ycoord=FEM.nodesCord(:,3);
Zcoord=FEM.nodesCord(:,4);
%%
%-----------------------------------------
Stru.length=max(abs(Xcoord));
Stru.width=max(abs(Ycoord));


%     Stru.thickness=1.04e-3;
%-----------D.O.F for each Node--------------
FEM.PlateNodeDof=5;%% for plate
Stiffener.nodedof=5; %% for stiffener
FEM.nodeCoordinates=FEM.nodesCord(:,2:3);

FEM.GDof=FEM.PlateNodeDof*size(FEM.nodeCoordinates,1);
FEM.typeplate='CQUAD8';
FEM.Dimension='2D';
switch FEM.typeplate
    case 'CQUAD4'
        FEM.GaussPointShear='1by1';
        FEM.GaussPointBend='2by2';
    case 'CQUAD8'
        FEM.GaussPointShear='2by2';
        FEM.GaussPointBend='3by3';
end
%---------------------------------------
FEM.NodeNumber=size(FEM.nodesCord,1);
FEM.elementNumber=size(FEM.elementNodes,1);


FEM.nodeCoordinates_label = zeros(size(FEM.nodeCoordinates,1),4);

FEM.nodeCoordinates_label(:,2:3) = FEM.nodeCoordinates;
figure(200);hold on;
plot(FEM.nodeCoordinates(:,1),FEM.nodeCoordinates(:,2),'ko')
FEM.nodeCoordinates_label(:,1)  = 1:size(FEM.nodeCoordinates,1);

%% plot
patch_plot(FEM.elementNodes,FEM.nodeCoordinates_label,200,'skin');axis image;


%% ================= Thermal TEMPERATURE  ================

Thermal.Temp_bot  = 1;
Thermal.k1 = 0.0;
Thermal.k2 = 0;
%% ================================================================
%% ===========-Composite Material Properties ======================
%% ================================================================
mat_type = 'GrEp';'CaEp';'GrEp';

Duran_Materials;

modeshape_fig_name ='TEST'  ;

%% 

Mat.kappa = 5/6;
Mat.G13=Mat.G12;
Mat.G23=Mat.G13;
Mat.v21=Mat.E2/Mat.E1*Mat.v12;
Mat.density=1800;

% to from a 3D orthotropic material
Mat.E3  = Mat.E2;
Mat.v13 = Mat.v12;
Mat.v31 = Mat.v13/Mat.E1*Mat.E3;
Mat.v23 = Mat.v12;
Mat.v32 = Mat.v23*Mat.E3/Mat.E2;

%
Laminate.layer= 4;
Laminate.layer_thickness = 1.016e-3/ Laminate.layer;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------PLATE STIFFNESS Matrix --------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FEM.typeplate='CQUAD8';

Laminate.thickness=Laminate.layer_thickness*ones(1,Laminate.layer);%% Uniform thickness
%
Stru.thickness= sum(Laminate.thickness);


center =[Stru.length Stru.width]/2;

physical_length = Stru.length; % change along width direction

% [MassPlate,MelemPlate] = LinearMassLaminatedPlate(FEM,Mat,Stru);

%% VAT
%           modeshape_fig_name ='buckling_modeshape_present_GraphteEpoxy'
t0  = T01(1);
t1  = T01(2);


flag = 'SYM';
half_number = Laminate.layer/2;
T0T1 = VAT_fiber(T01,half_number,flag);


[Kplate,Kelemp] = LinearStiffnessLaminatedPlate_VAT_v2_X_thermal(Mat, Stru,FEM,Laminate,T0T1,center(1),physical_length);


%% Thermal Load Analysis
FORCE_thermal = calculate_thermal_force(FEM,Laminate,Mat,Stru,T0T1,center,physical_length,Thermal);


%% -----------STATIC ANALYSIS ------------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORCE
% % distribute forces
%%Boundary Conditions
FEM.BCtype= 'SSSS-NXX-Thermal';%% 'SSSS-NYY'; 'SSSS-NXXNYY'; 'SSSS-NXY'

[ActiveDof,SidesDof]=InPlaneBCPlate5Dof(FEM);

FEM.displacement=zeros(1,FEM.GDof);
%                 FEM.displacement(ActiveDof)=Ktotal(ActiveDof,ActiveDof)\FORCE(ActiveDof);
FEM.displacement(ActiveDof)=Kplate(ActiveDof,ActiveDof)\FORCE_thermal(ActiveDof);
%
% % platetest; % static analysis of pure simply supported plate
ebar = 0;
postprocessISOTROPIC;

% static analysis for stress distributions

[FEM.stress,FEM.strain]=StressRecoveryPlate_VAT_center_average_v2_constant_angle(FEM,Laminate,Mat,Stru,T0T1,center,Stru.length);


[Thermal.stress,Thermal.strain]=StressRecoveryPlate_VAT_center_average_thermal_v1(FEM,Laminate,Mat,Stru,T0T1,center,physical_length,Thermal);


Thermal_Elastic.stress =FEM.stress - Thermal.stress;


FEM.stress =    Thermal_Elastic.stress;
% compute geometric stiffness due to thermal stress


KGplate =GeometryStiffnessPlate_stress_recovery(FEM,Stru,Laminate);


%% ============ Boundary Conditions ================
FEM.BCtype='SSSS-3'; %% Four simple supported sides

[ActiveDof,SideDof]=EssentialBCPlate5Dof(FEM);

%% =============Buckling Eigenvalue Analysis ================
Solver='buckling';
[V,D]=eigs(Kplate(ActiveDof,ActiveDof),-KGplate(ActiveDof,ActiveDof),10,'sm');
% % -----  buckling analysis of pure plate  ------

[DD,modeNo]=sort((diag(D)));
loadfactor=DD(1:10);
VVsort=V(:,modeNo);
frequency=loadfactor;
%-------------minimal load factor---------------
mineiglabel=find(min(abs(DD))==-DD);
if isempty(mineiglabel)==1
    mineiglabel=find(min(abs(DD))==DD);
end
plotmodenNo=mineiglabel;

%              modeshape_fig_name ='buckling_modeshape_present_GraphteEpoxy'

% plot mode shape
modeshapeplate;

critical_loadfactor = DD(mineiglabel);
DD(mineiglabel)


% ==============  END OF CODE  ============== %

% save('thermal_buckling_plate.mat')
