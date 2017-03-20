 function [kmean_out] = K_means_ding(I,K,centroids,space_dim,max_iter,region)
%% K means segmetation and SLIC segmentions (integrated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
%%%%% I         =  original image (RGB scale)
%%%%% K         =  original seeds number
%%%%% centroids =  user defined initial centers
%%%%% dimension =  3 refers to only RGB space, 5 refers x,y, RGB space
%%%%% maxiter   =  maximum iteration time, better to be less than 1000
%%%%% blk_len   =  in form the length of region S, if S = 0, process whole image
%%%%%              if S > 0 (in our case S = 50), The process region is 2Sx2S
%%%%%
%%%%%          Author: Houzhu Ding 04/02/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Startup
% Only for RGB image
if ~(size(I,3)==3)
    msg = 'Error input,only RGB.';
    error(msg);
end
I = double(I);
% Get image size
[m,n,channels]= size(I);
if space_dim == 5
    % Define the new mean position storage matrix
    x_mean_n      = zeros(1,K);
    y_mean_n      = zeros(1,K);
end
dist_wt = 4;%  = m^2 a factor that divide x,y by m=2 for SLIC
% Define the converge metric
converge      = 0;
% Define the cluster map, each pixel will be labeled 
clusters_map  = zeros(m,n)-1;
distance_map  = zeros(m,n)+inf;
% Define the final mean RGB storage matrix
R_mean_final  = zeros(1,K);
G_mean_final  = zeros(1,K);
B_mean_final  = zeros(1,K);
% Define the previous mean RGB storage matrix
R_mean_o      = zeros(1,K);
G_mean_o      = zeros(1,K);
B_mean_o      = zeros(1,K);
% Define the current mean RGB storage matrix
R_mean_n      = zeros(1,K);
G_mean_n      = zeros(1,K);
B_mean_n      = zeros(1,K);
% Define the distance between pixel to current mean
dist_p2om_set = zeros(1,K)+inf;

%% Step 1:Randomly pick K RGB triplets as seeds
% Define the current mean position storage matrix
x_mean_o       = centroids(1,:); 
y_mean_o       = centroids(2,:); 
% Calculate the mean RGB value of selected seeds
 for cnt = 1: K
     R_mean_o(cnt) = I(x_mean_o(cnt),y_mean_o(cnt),1);
     G_mean_o(cnt) = I(x_mean_o(cnt),y_mean_o(cnt),2);
     B_mean_o(cnt) = I(x_mean_o(cnt),y_mean_o(cnt),3);
 end
 %% Step 2:Given cluster centers, determine points in each cluster
 iter = 1;
 dist_om2nm_old = zeros(1,K)+inf;
 
 while ((converge) ~= 1 && iter < max_iter) 
     clusters_pt_num = zeros(1,K); % store pixel index of K clusters
     dist_om2nm     = zeros(1,K)+inf;
     fprintf('Iteration... %d\n',iter);
     % Count the cluster elements numbers
     % Calculate the new RGB sum value of each clauster
     R_sum = zeros(1,K);G_sum = zeros(1,K);B_sum = zeros(1,K);
     
     % Define the x,y space temp sum matrix to compute new x,y mean
     x_sum = zeros(1,K);y_sum = zeros(1,K);
     if region == 0
     %% Conventional K mean for Problem 1. Process region is whole image
         for i = 1:m
             for j = 1:n
                R = I(i,j,1);
                G = I(i,j,2);
                B = I(i,j,3);
                % Calculate K dist from one pixel to the current mean
                for cnt = 1: K 
                    % Default Space RGB for Problem 1
                     dist_p2m = sqrt((R-R_mean_o(cnt))^2+...
                                     (G-G_mean_o(cnt))^2+...
                                     (B-B_mean_o(cnt))^2);
                     dist_p2om_set(cnt) = dist_p2m; 
                     if space_dim == 5 % if x y should be considered recompute the distance
                     % Compute X Y RGB space distance for Problem 2 SLIC
                        dist_om2nm(cnt) =  sqrt(dist_om2nm(cnt)^2 +...
                                       (x_mean_n(cnt)-x_mean_o(cnt))^2/dist_wt+...
                                       (y_mean_n(cnt)-y_mean_o(cnt))^2/dist_wt);
                     end
                end
                % Sort the dist and the minimum dist order is its cluster
                [min_dist,min_order] = min(dist_p2om_set);
                clusters_map(i,j) = min_order;
             end
         end
     else
     %% SLIC: Optional compare the piexl within 100(2xblock size) pixels
         for cnt = 1:K
             for i = ceil(x_mean_o(cnt)-region):ceil(x_mean_o(cnt)+region)
                 for j = ceil(y_mean_o(cnt)-region):ceil(y_mean_o(cnt)+region)                
                     if i > 0 && j > 0 && i < m && j < n
                        R = I(i,j,1);
                        G = I(i,j,2);
                        B = I(i,j,3);
                        dist_p2m = sqrt((R-R_mean_o(cnt))^2+...
                                 (G-G_mean_o(cnt))^2+...
                                 (B-B_mean_o(cnt))^2);
                         % Compute distance in 5D space
                        if space_dim == 5
                            dist_p2m = sqrt(dist_p2m^2+...
                                        (i - (x_mean_o(cnt)))^2/dist_wt+...
                                        (j - (y_mean_o(cnt)))^2/dist_wt); 
                        end
                        if dist_p2m < distance_map(i,j)
                            distance_map(i,j) = dist_p2m;
                            clusters_map(i,j) = cnt;
                        end
                     end
                 end
             end
         end
     end  
     % Compute the summary of each channel or position
     for i = 1:m
         for j = 1:n   
             for cnt = 1: K 
                if clusters_map(i,j) == cnt
                    clusters_pt_num(cnt) = clusters_pt_num(cnt) + 1; 
                    % Add all elements in each RGB channel/position
                    R_sum(cnt) = R_sum(cnt) + I(i,j,1);
                    G_sum(cnt) = G_sum(cnt) + I(i,j,2);
                    B_sum(cnt) = B_sum(cnt) + I(i,j,3);
                    x_sum(cnt) = x_sum(cnt) + i;
                    y_sum(cnt) = y_sum(cnt) + j;
                end
             end
         end
     end
%% Calculate the new cluster centers and update
    % new mean and previous mean
    for cnt = 1: K
        if clusters_pt_num(cnt)
            R_mean_n(cnt) = R_sum(cnt)/clusters_pt_num(cnt);
            G_mean_n(cnt) = G_sum(cnt)/clusters_pt_num(cnt);
            B_mean_n(cnt) = B_sum(cnt)/clusters_pt_num(cnt);
            % Compute new centroids
            x_mean_n(cnt) = x_sum(cnt)/clusters_pt_num(cnt);
            y_mean_n(cnt) = y_sum(cnt)/clusters_pt_num(cnt);
            % Compute RGB space distance
            dist_om2nm(cnt) = sqrt((R_mean_n(cnt)-R_mean_o(cnt))^2+...
                              (G_mean_n(cnt)-G_mean_o(cnt))^2+...
                              (B_mean_n(cnt)-B_mean_o(cnt))^2);   
            if space_dim == 5 % if x y should be considered recompute the distance
                % Compute X Y RGB space distance
                dist_om2nm(cnt) =  sqrt(dist_om2nm(cnt)^2 +...
                                   (x_mean_n(cnt)-x_mean_o(cnt))^2/dist_wt+...
                                   (y_mean_n(cnt)-y_mean_o(cnt))^2/dist_wt);
            end
        end
    end
    % Coverage condition: the cluster centers don't change
    if (dist_om2nm_old == dist_om2nm)
        converge = 1;
        % Record the mean RGB value 
        R_mean_final = R_mean_n;
        G_mean_final = G_mean_n;
        B_mean_final = B_mean_n;
    else
            % Update the new mean pixels
            dist_om2nm_old = dist_om2nm;
            R_mean_o = R_mean_n;
            G_mean_o = G_mean_n;
            B_mean_o = B_mean_n;
            x_mean_o = x_mean_n;
            y_mean_o = y_mean_n;
    end
    % Count the iteration times no more than 100 (usr defined)
    iter = iter + 1;
 end
%% Step 3 : Label clusters with average RGB
for i = 1:m
    for j = 1:n
        for cnt = 1: K
            if (clusters_map(i,j) == cnt)
                I(i,j,:) = [R_mean_final(cnt),G_mean_final(cnt),B_mean_final(cnt)];
            end
        end
    end
end
% Show black boundary for problem 2 SLIC algorithm
if (space_dim == 5)              
    for i = 2:m
        for j = 2:n
            if clusters_map(i,j) ~= clusters_map(i-1,j)
                I(i,j,:) = [0,0,0];
            else
                if  clusters_map(i,j) ~= clusters_map(i,j-1)
                    I(i,j,:) = [0,0,0];
                end
            end
            
        end
    end
    for i = 1
        for j = 2:n
            if  clusters_map(i,j) ~= clusters_map(i,j-1)
                I(i,j,:) = [0,0,0];
            end
        end
    end
    
    for i = 2:m
        for j = 1
            if  clusters_map(i,j) ~= clusters_map(i-1,j)
                I(i,j,:) = [0,0,0];
            end
        end
    end
end
% Convert the image back to unit8 format
kmean_out = uint8(I);
% imshow(uint8(I));
%% The end of K - mean algorithm /SLIC algorithm
     