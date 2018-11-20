%% Material Constant

function [bendmatrix,shearmatrix, stretchmatrix]=materialplate(Mat,Stru)
    display('Isotropic Material-Aluminum')
    E=Mat.E;possion=Mat.possion;G=Mat.G; kappa=Mat.kappa;
    % Material Constant D matrix \\ {sigma}=[D]{epsilon};
    % Uniform stiffness
    I=1/12*Stru.thickness.^3;
    % constitive law matrix [D]
    D_bend=E/(1-possion^2)*[1,possion,0;possion,1,0;0,0,(1-possion)/2];
    D_shear=[G,0;0,G];
    % [D]--constitutive law matrix constant
    D=[E/(1-possion^2),E*possion/(1-possion^2),0,0,0;...
        E*possion/(1-possion^2),E/(1-possion^2),0,0,0;...
        0,0,E/(1-possion^2)*(1-possion)/2,0,0;...
        0,0,0,G,0;...
        0,0,0,0,G];
    % C_Bending; C_Shearing; Bending and shearing stiffness;
    bendmatrix=I*D_bend;
    shearmatrix=kappa*Stru.thickness*D_shear;
    stretchmatrix=Stru.thickness*D_bend;