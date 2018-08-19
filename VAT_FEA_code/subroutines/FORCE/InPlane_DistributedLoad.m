function p=InPlane_DistributedLoad(x,type,alpha,b)

% x     ---- location of nodes where the load applied;
% type  ---- linear or nonlinear load distribution;
% alpha ---- for different linear variations;
% b     ---- length of applied load side.

switch type 
  
    case 'linear'
        
        p=1-alpha*x/b;

    case 'quadratic'
        
        p=4*x/b-4*x^2/b^2;
       
end
