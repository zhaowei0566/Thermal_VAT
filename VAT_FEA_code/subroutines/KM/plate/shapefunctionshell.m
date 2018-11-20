function [shape,naturalderivatives,d2Nds2]=shapefunctionshell(xi,eta,FEM)
type=FEM.typeplate;
switch type
    
    case 'CQUAD4'
        shape=1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);(1+xi)*(1+eta);(1-xi)*(1+eta)];
        % naturalderivatives=partial shape w.r.t xi or eta
        
        naturalderivatives=1/4*[-(1-eta),(1-eta),(1+eta),-(1+eta);...
            -(1-xi),-(1+xi),(1+xi),(1-xi)];
        
    case 'CQUAD8'
        shape=[1/4*(1-xi)*(1-eta)*(-xi-eta-1);
            1/4*(1+xi)*(1-eta)*(xi-eta-1);
            1/4*(1+xi)*(1+eta)*(xi+eta-1);
            1/4*(1-xi)*(1+eta)*(-xi+eta-1);
            1/2*(1-eta)*(1+xi)*(1-xi);
            1/2*(1+xi)*(1+eta)*(1-eta);
            1/2*(1+eta)*(1+xi)*(1-xi);
            1/2*(1-xi)*(1+eta)*(1-eta)];
        naturalderivatives(1,:)=[ - (xi/4 - 1/4)*(eta - 1) - ((eta - 1)*(eta + xi + 1))/4
            ((eta - 1)*(eta - xi + 1))/4 - (xi/4 + 1/4)*(eta - 1)
            (xi/4 + 1/4)*(eta + 1) + ((eta + 1)*(eta + xi - 1))/4
            (xi/4 - 1/4)*(eta + 1) + ((eta + 1)*(xi - eta + 1))/4
            (eta/2 - 1/2)*(xi - 1) + (eta/2 - 1/2)*(xi + 1)
            -((eta - 1)*(eta + 1))/2
            - (eta/2 + 1/2)*(xi - 1) - (eta/2 + 1/2)*(xi + 1)
            ((eta - 1)*(eta + 1))/2];
        naturalderivatives(2,:)=[ - (xi/4 - 1/4)*(eta - 1) - (xi/4 - 1/4)*(eta + xi + 1)
            (xi/4 + 1/4)*(eta - xi + 1) + (xi/4 + 1/4)*(eta - 1)
            (xi/4 + 1/4)*(eta + 1) + (xi/4 + 1/4)*(eta + xi - 1)
            (xi/4 - 1/4)*(xi - eta + 1) - (xi/4 - 1/4)*(eta + 1)
            ((xi - 1)*(xi + 1))/2
            - (xi/2 + 1/2)*(eta - 1) - (xi/2 + 1/2)*(eta + 1)
            -((xi - 1)*(xi + 1))/2
            (xi/2 - 1/2)*(eta - 1) + (xi/2 - 1/2)*(eta + 1)];

end
