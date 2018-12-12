function [contX, contY, contAreas] = MedianArea(Im);
%{
    Get median area for each candidate thresh
    Args:
        Im - Gray-scale image

    Returns:
        contX - X-values of contours for each identified nuclei
        contY - Y-values of contours for each identified nuclei
        contAreas - Areas for each contour
    
    Ref: https://stackoverflow.com/questions/20040661/finding-2d-area-defined-by-contour-lines-in-matlab
%}
%figure, 
[C, h] = imcontour(Im, 1); 

if size(C,2) == 0 % No contours!
    contAreas = 0;
    contX = 0;
    contY = 0;
    return;
end

C_x = C(1,:);
C_y = C(2,:);
% Find start pts for each contour
start(1) = 1;
k = 1;
while start(k) <= length(C_x)
    new_start = start(k) + C_y(start(k)) + 1; 
    if new_start > length(C_x)
        break
    end
    start(k+1) = new_start; 
    k = k+1;
end

% Calc area
for i = 1:length(start)
    len = C_y(start(i));
    contX{i} = C_x(start(i)+1:start(i)+len);
    contY{i} = C_y(start(i)+1:start(i)+len);
    contAreas(i) = polyarea(contX{i}, contY{i});
end