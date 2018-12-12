imagePath = 'test.tif'; % path + name of file
imageName = 'test_enh'; % output file name
A = imread(imagePath);

% Contrast enhancement
R = A(:,:,1);
adj = imadjust(R);

%% Object-based detection
gray = adj;
% Get PWM curve
pwm = PWM(gray);
% Fit 15 poly function
pwm_x = 0:255;
p = polyfit(pwm_x, pwm, 15);
p_x = 0:255;
p_y = polyval(p, p_x);

% -- Find candidate thresholds 
% Create poly function
syms x
f = p(1)*x^(15) + p(2)*x^(14) + p(3)*x^(13) + p(4)*x^(12) + p(5)*x^(11) ...
  + p(6)*x^(10) + p(7)*x^(9)  + p(8)*x^(8)  + p(9)*x^(7)  + p(10)*x^(6) ...
  + p(11)*x^(5) + p(12)*x^(4) + p(13)*x^(3) + p(14)*x^(2) + p(15)*x^(1) + p(16);
f1 = diff(f); % 1st deriv
f2 = diff(f1); % 2nd deriv
% Set 2nd deriv equal to 0
inflec_pt = solve(f2, 'MaxDegree', 15);
inflec_pt = double(inflec_pt);
inflec_pt(inflec_pt > 255) = []; % remove nums out of bounds

% Get initial threshold
maxMedArea = 0;
for i = 1:length(inflec_pt)
%    i
    temp = gray;
    if imag(inflec_pt(i)) ~= 0
	continue
    end
    temp(gray <= inflec_pt(i)) = 0;
    temp(gray > inflec_pt(i)) = 255;
    [contX, contY, contAreas] = MedianArea(temp);
    if contAreas == 0 % No contours!
        continue;
    end
    tempMedArea = median(contAreas);
    if tempMedArea > maxMedArea
        maxMedArea = tempMedArea;
        maxContAreas = contAreas;
        maxContX = contX;
        maxContY = contY;
        maxThresh = inflec_pt(i);
    end
end

%% Area-based correction
% Set new variables, used for testing purposes
Im = gray;
contX = maxContX;
contY = maxContY;
initialThresh = maxThresh;
contAreas = maxContAreas;

% Identify small, normal, and big nuclei
meanArea = mean(contAreas);
smallThresh = 0.25*meanArea;
bigThresh = meanArea*5;

tooSmall = [];
tooBig = [];
for i = 1:length(contAreas) % Get indices of nuclei that are too small or big
    if contAreas(i) < smallThresh
        tooSmall = [tooSmall i];
    elseif contAreas(i) > bigThresh
        tooBig = [tooBig i];
    end
end

% Resolve big nuclei, decrease initial thresh
for i = 1:length(tooBig) % get contours for big nuclei 
    idx = tooBig(i);
    bigContoursX{i} = maxContX{idx};
    bigContoursY{i} = maxContY{idx};
end

% Change initial threshold iteratively
allBigX = {};
allBigY = {};
stepSize = 2;
for i = 1:length(bigContoursX) 
    % Get mask for nuclei cluster
    bigMask = poly2mask(bigContoursX{i}, bigContoursY{i}, size(Im,1), size(Im,2));
    bigMask = uint8(bigMask);
    % Apply mask
    J = Im.*bigMask;
    % Convert background (complement of objects found) to white (255)
    for i = 1:size(Im,1)
        for j = 1:size(Im,2)
            if bigMask(i,j) == 0
                J(i,j) = 255;
            end
        end	
    end

    for j = 1:stepSize:initialThresh 
        clear contX contY contAreas
        i, j
        newThresh = initialThresh - j;
        % Binarize w/ new threshold
        temp = J;
        temp(J <= newThresh) = 0;
        temp(J > newThresh) = 255;
        % Get new contours
        [contX, contY, contAreas] = MedianArea(temp);
        % Check if all nuclei areas are no longer big
        if all(contAreas < bigThresh) == 1 && all(contAreas > smallThresh) == 1
            break;
        end
    end
    % Record all contours
    allBigX = [allBigX, contX]; 
    allBigY = [allBigY, contY];
end

% Update mask, take out big and small nuclei 
restContX = {};
iter = 1;
for i = 1:length(maxContX)
    if ismember(i, tooBig) == 1 
        continue;
    end
    restContX{iter} = maxContX{i};
    restContY{iter} = maxContY{i};
    iter = iter + 1;
end

% Combine all contours = normal sized nuclei + big nuclei broken down
combContX = [allBigX, restContX];
combContY = [allBigY, restContY];

% Make combined mask
totMask = zeros(size(Im));
for i = 1:length(combContX)
    nextMask = poly2mask(combContX{i}, combContY{i}, size(Im,1), size(Im,2));
    totMask = nextMask + totMask; 
end

%% Nuclei Separation

% Marker-controlled watershed algorithm
% https://journals.sagepub.com/doi/full/10.1155/2013/268046
dist = bwdist(~totMask); % Distance transform
% Normalize
dist_norm = dist - min(dist(:));
dist_norm = dist_norm / max(dist_norm(:));

ext = imextendedmax(dist_norm, 0.2); % Extended-max transform

% Watershed
ext2 = ext;
ext2(~dist_norm) = Inf;
watersh = watershed(ext2);
watersh(watersh~=0) = 1;
 
% Apply watershed to mask
totMask_ws = uint8(watersh) .* uint8(totMask);
totMask_ws = ~(totMask .* 255);

% Get contours for new mask
[wsContX, wsContY, wsContAreas] = MedianArea(totMask_ws);

wsContX_filt = wsContX;
wsContY_filt = wsContY;

% Get center coordinates of found nuclei
for i = 1:length(wsContX_filt)
    centerX(i) = mean(wsContX_filt{i});
    centerY(i) = mean(wsContY_filt{i});
end

% Write centers to file
fileID = fopen(strcat('centers_',imageName,'.txt'), 'w');
fprintf(fileID, '%f %f\n', [size(A,1); size(A,2)]);
fprintf(fileID, '%f %f\n', [centerX; centerY]);
