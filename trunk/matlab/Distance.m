function [ d ] = Distance( p1, p2 )
%Distance assess the euclidean distance between p1 and p2 points
%   p1 = [x,y], p2 = [x,y]
    d = sqrt( (p2(1)-p1(1))^2 + (p2(2)-p1(2))^2 );
end

