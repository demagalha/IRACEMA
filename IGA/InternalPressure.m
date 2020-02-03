function force = InternalPressure(x,y,z,R)
if (x^2 +y^2) < R+eps
    theta = atan(y/x);
    a = cos(theta);
    b = sin(theta);
    c = 0;
else
    a = 0; b = 0; c = 0;
end
force = [a; b; c];
end