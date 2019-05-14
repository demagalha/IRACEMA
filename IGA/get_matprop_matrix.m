function D = get_matprop_matrix(idx, YOUNG, POISSON)
% 1 = Isotropic
% 2 = Ortotropic
% 3 = Anisotropic
switch idx
    case 1
        lambda = YOUNG*POISSON/((1+POISSON)*(1-2*POISSON)); % Lamé
        mi = YOUNG/(2*(1+POISSON)); % G
        D = [lambda+2*mi, lambda,      lambda,      0,   0,   0;
             lambda,     lambda+2*mi,  lambda,      0,   0,   0;
             lambda,     lambda,      lambda+2*mi,  0,   0,   0;
             0,          0,           0,           mi,  0,   0;
             0,          0,           0,           0,   mi,  0;
             0,          0,           0,           0,   0,   mi;];
    case 2
        D = [];
    case 3
        D = [];
            
end