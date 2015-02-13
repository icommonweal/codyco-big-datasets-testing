function [ ecc ] = ellipsoid_ecc( p )
%ELLIPSOID_ECC Compute the two eccentrisies of an ellipoisd
%   Takes in input the ellipsoid in implicit form
    [c,ax] = ellipsoid_im2ex(p);
%     ecc(1) = ax(2)/ax(1);
%     ecc(2) = ax(3)/ax(1);
%     ecc(3) = sqrt(1-(ax(3)^2)/(ax(1)^2));   
%     ecc(4) = sqrt(1-(ax(3)^2)/(ax(2)^2));
    ecc = ax;
end

