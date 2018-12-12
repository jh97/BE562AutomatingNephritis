function [ area ] = findarea( h )
%% Finds the area enclosed by a contour line
   % h - The Contour object created using the contour() function
 	xx = h.ContourMatrix(1,:);
 	yy = h.ContourMatrix(2,:);
 	area(1) = polyarea(xx,yy);
end