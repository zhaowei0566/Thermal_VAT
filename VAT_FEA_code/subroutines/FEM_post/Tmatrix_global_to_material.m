function T = Tmatrix_global_to_material(theta)


m = cosd(theta);
n = sind(theta);

T = [m^2 n^2 2*m*n;
    n^2 m^2 -2*m*n;
    -m*n m*n m^2-n^2];