%% Material Constant
function [Amatrix,Dmatrix,Ashear,Bmatrix,Qbar_theta,ThermalExpCoeff]=LaminatedComposite(Mat,Stru,Laminate,layerNo)
format long;
theta=Laminate.theta(layerNo);
thickness=Laminate.thickness(layerNo);
display(['Laminated Composites Plate - Orthotropic Materials, Layer=' num2str(layerNo)])
c=cos(theta/180*pi);
s=sin(theta/180*pi);
%--------Transform Matrix for Coordinate Transformation--------
TransformMatrix=[c^2,s^2,0,0,0,2*s*c;...
    s^2,c^2,0,0,0,-2*s*c;...
    0,0,1,0,0,0;...
    0,0,0,c,-s,0;...
    0,0,0,s,c,0;...
    -s*c,s*c,0,0,0,c^2-s^2]; %% Global to Principal of each layer


Q11=Mat.E1/(1-Mat.v12*Mat.v21);
Q12=Mat.v12*Mat.E2/(1-Mat.v12*Mat.v21);
Q22=Mat.E2/(1-Mat.v12*Mat.v21);
Q66=Mat.G12;
Q44=Mat.G23;
Q55=Mat.G13;

Q_inplane=[Q11 Q12 0   0   0;
           Q12 Q22 0   0   0;
            0   0 Q44  0   0;
            0   0  0  Q55  0;
            0   0  0   0  Q66];

%------Orthotropic Material--------------
Q11b=Q11*c^4+2*(Q12+2*Q66)*s^2*c^2+Q22*s^4;
Q12b=(Q11+Q22-4*Q66)*s^2*c^2+Q12*(s^4+c^4);
Q22b=Q11*s^4+2*(Q12+2*Q66)*s^2*c^2+Q22*c^4;
Q16b=(Q11-Q12-2*Q66)*s*c^3+(Q12-Q22+2*Q66)*s^3*c;
Q26b=(Q11-Q12-2*Q66)*c*s^3+(Q12-Q22+2*Q66)*c^3*s;
Q66b=(Q11+Q22-2*Q12-2*Q66)*s^2*c^2+Q66*(s^4+c^4);
Q44b=Q44*c^2+Q55*s^2;
Q45b=(Q55-Q44)*c*s;
Q55b=Q55*c^2+Q44*s^2;


Qbar=[Q11b,  Q12b,  0,   0,  Q16b;
      Q12b,  Q22b,  0,   0,  Q26b;
       0,     0,  Q44b, Q45b,  0;
       0,     0,  Q45b, Q55b,  0;
     Q16b,   Q26b, 0,    0,   Q66b];
 
bendingNO=[1,2,5];shearNO=[3,4];
number=[1,2,3,4,5];
Qbar_theta=Qbar(number,number);

D_bend=Qbar(bendingNO,bendingNO);
D_shear=Qbar(shearNO,shearNO);

zk1=-Stru.thickness/2+layerNo*thickness;
zk=-Stru.thickness/2+(layerNo-1)*thickness;
Amatrix=(zk1-zk)*D_bend;
Bmatrix=1/2*(zk1^2-zk^2)*D_bend;
Dmatrix=1/3*(zk1^3-zk^3)*D_bend;
Ashear=Mat.kappa*(zk1-zk)*D_shear;

alphaX=Mat.alpha1*c^2+Mat.alpha2*s^2;
alphaY=Mat.alpha1*s^2+Mat.alpha2*c^2;
alphaXY=2*(Mat.alpha1-Mat.alpha2)*c*s;
 
ThermalExpCoeff=[alphaX alphaY alphaXY]; 

%% End of diary