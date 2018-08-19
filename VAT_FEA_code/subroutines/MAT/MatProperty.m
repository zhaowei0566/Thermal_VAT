%% Chattopadhady 1 curvilinear stiffener

%     Mat.kappa=5/6;
%
%     Mat.E1=19.2e6*psipa;
%     Mat.E2=1.56e6*psipa;
%
%     Mat.G12=0.82e6*psipa;
%     Mat.G13=Mat.G12;
%     Mat.G23=0.49e6*psipa;
%
%     Mat.v12=0.24;
%     Mat.v21=Mat.E2/Mat.E1*Mat.v12;
%
%     Mat.density=1800;
%
%     Mat.alpha1=-0.04e-6;
%     Mat.alpha2=16.7e-6;

%% Material in Mitt's paper

% 
% Mat.kappa=5/6; %% Mitt..
% Mat.E1=138e9;
% Mat.E2=8.96e9;
% % Mat.G12=7.372e9;
% Mat.G12=7.1e9;
% Mat. G13=Mat.G12;Mat.G23=Mat.G12;
% % Mat.G12=7.372e9; Mat. G13=Mat.G12;Mat.G23=Mat.G12; % for shear load case
% Mat.v12=0.30;
% Mat.v21=Mat.E2/Mat.E1*Mat.v12;
% Mat.density=1800;
% Mat.alpha1=-0.04e-6;
% Mat.alpha2=16.7e-6;

%% Houdayfa Thermal Example
% Mat.E1=22.5e6;
% Mat.E2=1.17e6;
% Mat.G12=0.66e6;Mat.G13=Mat.G12;
% Mat.G23=Mat.G12;
% Mat.v12=0.22;
% Mat.alpha1=-0.04e-6;
% Mat.alpha2=16.7e-6;
% Mat.density=1800;

%%  PARAMETRIC STUDIES MATERIAL PROPERTIES %%%
%     Mat.possion=0.3; %% For Stiffener dimension
%     Mat.density=2823;
%     RR=40; %%orthogonic ratio
%     Mat.kappa=5/6;
%     Mat.E2=1e7;
%     Mat.E1=RR*Mat.E2;
%
%     Mat.G12=0.6*Mat.E2; Mat. G13=Mat.G12;
%     Mat.G23=0.5*Mat.E2;
%     Mat.v12=0.25;

%% Lee's vibration case material

Mat.kappa=5/6;

Mat.E1=128e9;
Mat.E2=11e9;

Mat.G12=4.48e9;
Mat.G13=4.48e9;
Mat.G23=1.53e9;

Mat.v12=0.25;
Mat.v21=Mat.E2/Mat.E1*Mat.v12;

Mat.density=1500;

Mat.alpha1=-0.04e-6;
Mat.alpha2=16.7e-6;


