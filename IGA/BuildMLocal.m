function Me = BuildMLocal(R,Jmod,rho)
N = zeros(length(R),3);
    for i=1:length(R)
        N(1:3,1+3*(i-1):3+3*(i-1)) = R(i)*eye(3);
    end
    Me = Jmod*rho*N'*N;
end