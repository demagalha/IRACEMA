function force = PinchedForce(x,y,z,R)
        force = [0; 0; 0];
    if (z >= 0) && (z <= 3)
        if (x >= R-0.01*R) && (x <= R)
            force = [-1; 0; 0];
        end
    end
end