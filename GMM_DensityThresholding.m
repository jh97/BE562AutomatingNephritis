test_image = 'test.tif';
test_coord = 'centers_test.txt';
threshold_density = 0.00639;
mean_thresh = 0.05;

fileID = fopen(test_coord,'r'); 
sizeX = [2 Inf];
formatSpec = '%f %f';

X = fscanf(fileID,formatSpec,sizeX);
X = X';
image_dim = X(1,:);
X = X(2:end,:);
image_area = image_dim(1)*image_dim(2);
%% STEP 1: Choose initial values for our parameters.

N = size(X, 1); % Number of data points
K = round(image_area/250000, 0);  % Estimate number of clusters to generate.

% Randomly select k data points to serve as the initial means.
indeces = randperm(N);
mu = X(indeces(1:K), :);

cov_matrix = [];

% Use the overall covariance of the dataset as the initial variance for each cluster.
for (k = 1 : K)
    cov_matrix{k} = cov(X);
end

% Assign equal prior probabilities to each cluster.
alpha = ones(1, K) * (1 / K);

%%===================================================
%% STEP 2: Run Expectation Maximization

% Matrix to hold the probability that each data point belongs to each cluster.
% One row per data point, one column per cluster.
w = zeros(N, K);

% Loop until convergence.
for (i = 1:1000)
    
    fprintf('  EM Iteration %d\n', i);
    %% Expectation
    %
    % Calculate the probability for each data point for each distribution.
    
    % Matrix to hold the pdf value for each every data point for every cluster.
    % One row per data point, one column per cluster.
    pdf = zeros(N, K);
    
    % For each cluster...
    for k = 1 : K       
        % Evaluate the bivariate gaussian pdf for all data points for cluster 'k'.
        pdf(:, k) = gaussianND(X, mu(k, :), cov_matrix{k});
    end
    
    % Multiply each pdf value by the prior probability for cluster.
    pdf_w = bsxfun(@times, pdf, alpha);
    
    % Divide the weighted probabilities by the sum of weighted probabilities across each cluster.
    w = bsxfun(@rdivide, pdf_w, sum(pdf_w, 2));
    
    %% Maximization
    % Store the previous means.
    prevMu = mu;    
    
    % For each of the clusters...
    for k = 1 : K
    
        % Calculate the prior probability for cluster 'k'.
        alpha(k) = mean(w(:, k), 1);
        
        % Calculate the new mean for cluster 'k' by taking the weighted
        % average of all data points.
        mu(k, :) = weightedAverage(w(:, k), X);

        % Calculate the covariance matrix for cluster 'k' by taking the 
        % weighted average of the covariance for each training example. 
        cov_k = zeros(2, 2);
        
        % Subtract the cluster mean from all data points.
        Xm = bsxfun(@minus, X, mu(k, :));
        
        % Calculate the contribution of each training example to the covariance matrix.
        for i = 1 : N
            cov_k = cov_k + (w(i, k) .* (Xm(i, :)' * Xm(i, :)));
        end
        
        % Divide by the sum of weights within cluster 'k'.
        cov_matrix{k} = cov_k ./ sum(w(:, k));
    end
    
    % Check for convergence.
    if (abs(mu-prevMu) <= mean_thresh)
        break
    end  
end

%% STEP 3: Plot the data points and their estimated contour lines.
figure(1);
imshow(test_image);
daspect([1 1 1]);
hold on;
plot(X(:, 1), X(:, 2), 'g.', 'MarkerSize', 5);

set(gcf,'color','white') % White background for the figure.

 % Create a [1,000,000 x 2] matrix 'grid' of coordinates representing
 % the input values over the grid.
 gridSize = 1000;
 points = linspace(0, max(image_dim(1,:)), gridSize);
 [A B] = meshgrid(points, points);
 grid = [A(:), B(:)];
 
 density = [];
  
 for i = 0:(K-1)
 % Calculate the Gaussian response for every value in the grid.
 z = gaussianND(grid, mu(i+1, :), cov_matrix{i+1});

 % Reshape the responses back into a 2D grid to be plotted with contours.
 Z(i*gridSize+1:(i+1)*gridSize, 1:gridSize) = reshape(z, gridSize, gridSize);
  
 circles = 1; % Number of contours to plot
 % Plot a contour line for each cluster to show its pdf over the data.
 [C1, h1] = contour(points, points, Z(i*gridSize+1:(i+1)*gridSize, 1:gridSize), circles, 'LineWidth', 2, 'LineColor', 'r');
 C1 = C1(:,2:end);
 if (isempty(C1))
    density(i+1) = 0;
 else %% Calculate the cell density within each cluster
    area1 = findarea(h1);
    [in1, on1] = inpolygon(X(:,1),X(:,2),C1(1,:),C1(2,:));
    total1 = sum(in1) + sum(on1); 
    density(i+1) = total1/area1;
 end
 end
 
 axis ij
 axis([0 image_dim(1,2) 0 image_dim(1,1)])
 set(gca,'DataAspectRatio',[1 1 1])
 
 %% STEP 4: Perform Density Thresholding
 figure(2);
 imshow(test_image)
 daspect([1 1 1]);
 hold on;
 plot(X(:, 1), X(:, 2), 'g.', 'MarkerSize', 5);
 set(gcf,'color','white') % White background for the figure.

 % Only plot contours with regions of nephritis
 for i = 0:(K-1)
 if (density(i+1) >= threshold_density)
     contour(points, points, Z(i*gridSize+1:(i+1)*gridSize, 1:gridSize), circles, 'LineWidth', 2, 'LineColor', 'r');
 end
 end
 
 axis ij
 axis([0 image_dim(1,2) 0 image_dim(1,1)])
 set(gca,'DataAspectRatio',[1 1 1])