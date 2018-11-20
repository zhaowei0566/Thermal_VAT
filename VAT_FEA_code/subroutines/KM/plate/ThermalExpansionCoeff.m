function [alphaX,alphaY,alphaXY] = ThermalExpansionCoeff(Mat)

c=cos(theta/180*pi);
s=sin(theta/180*pi);


alphaX = Mat.alpha1*c^2 + Mat.alpha2*s^2;


alphaY = Mat.alpha1*s^2 + Mat.alpha2*c^2;



alphaXY = 2 *(Mat.alpha1- Mat.alpha2)*s*c;

