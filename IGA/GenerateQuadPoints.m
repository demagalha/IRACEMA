function [quad_point_index, weights] = GenerateQuadPoints(poly_orders)
    switch length(p)
        case 1
            [quad_point_index, weights] = getGP(p(1));
        case 2
            [qu, wu] = getGP(p(1));
            [qv, wv] = getGP(p(2));
            N_QUAD = zeros(2,1);
            N_QUAD(1) = length(qu);
            N_QUAD(2) = length(qv);
            N_QUAD = prod(N_QUAD);
            quad_point_index = zeros(N_QUAD,2);
            weights = zeros(N_QUAD,1);
            for n = 1:N_QUAD
                [i,j] = ind2sub([length(qu),length(qv)],n);
                quad_point_index(n,:) = [qu(i),qv(j)];
                weights(n) = wu(i)*wv(j);
            end
        case 3
            [qu, wu] = getGP(p(1));
            [qv, wv] = getGP(p(2));        
            [qw, ww] = getGP(p(3));        
            N_QUAD = zeros(3,1);
            N_QUAD(1) = length(qu);
            N_QUAD(2) = length(qv);
            N_QUAD(3) = length(qw);
            N_QUAD = prod(N_QUAD);
            quad_point_index = zeros(N_QUAD,3);
            weights = zeros(N_QUAD,1);
            for n = 1:N_QUAD
                [i,j,k] = ind2sub([length(qu),length(qv),length(qw)],n);
                quad_point_index(n,:) = [qu(i),qv(j),qw(k)];
                weights(n) = wu(i)*wv(j)*ww(k);
            end
    end
end