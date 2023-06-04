%%%%%%%%%%%%%%%%%%%%% Part 1&2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
tic
%% Read Input Image and convert to CieLAB Space

[file,path] = uigetfile('*.*');
f = fullfile(path,file);
a = imread(f);
% a = imresize(a,0.5);
a_lab = rgb2lab(a);
%% Uncomment for the fire effect
% labTransformation = makecform('srgb2lab');
% a_lab = double(applycform(a,labTransformation));

%% Parameters
% 在這裡嘗試不同m,k來觀察結果的變化
m = 20;
n = 1;
k = 10;
%%

m = 15;
n = 5;
k = 200;
%% 
N = size(a,1)*size(a,2);

s = sqrt(N/k);

%% Gradient Image
G = zeros(size(a,1)-1,size(a,2)-1);
for i = 2:size(a,1)-1
    for j = 2:size(a,2)-1
        gx = (squeeze(a_lab(i+1,j,:))-squeeze(a_lab(i-1,j,:)));
        gy = (squeeze(a_lab(i,j+1,:))-squeeze(a_lab(i,j-1,:)));
        G(i,j) = gx(1)^2 + gx(2)^2 + gx(3)^2 + gy(1)^2 + gy(2)^2 + gy(3)^2;
    end
end
% figure;
% imagesc(G);

%% Initializing the Centers
s = ceil(s);
cx = s:s:size(a,1)-s;
cy = s:s:size(a,2)-s;
p=1;
for i = 1:size(cx,2)
    for j = 1:size(cy,2)
        loc(p,:) = [cx(i),cy(j)];
        p=p+1;
    end
end

for i = 1:size(loc,1)
    c(i,:) = [a_lab(loc(i,1),loc(i,2),1) a_lab(loc(i,1),loc(i,2),2) a_lab(loc(i,1),loc(i,2),3) loc(i,1) loc(i,2)];
end

% %% SLIC Algorithm (uncomment for correct implementation of the SLIC
% algortihm. the current setup emulates what would happen if this
% initialization step is skipped.
% win = 7;
% n1 = floor(win/2);
% 
% lochange = -n1:n1;
% 
% 
% for i = 1:size(loc,1)
%     H = G(loc(i,1)-n1:loc(i,1)+n1,loc(i,2)-n1:loc(i,2)+n1);
%     [a1,b1] = min(H);
%     [a2,b2] = min(a1);
%     loc(i,1) = loc(i,1) + lochange(b1(b2));
%     loc(i,2) = loc(i,2) + lochange(b2);
%     c(i,:) = [a_lab(loc(i,1),loc(i,2),1) a_lab(loc(i,1),loc(i,2),2) a_lab(loc(i,1),loc(i,2),3) loc(i,1) loc(i,2)];
% end

iter = 0;
msg = 'Segmenting ...';
x = 0;
f = waitbar(x,msg);
while iter < n
   
   for i2 = 1:size(a,1)
       for j2 = 1:size(a,2)
           dis = [];
           for k2 = 1:size(loc,1)
               if sqrt((i2-loc(k2,1))^2 + (j2 - loc(k2,2))^2) < 2*s
                   d = sqrt((a_lab(i2,j2,1)-c(k2,1))^2 + (a_lab(i2,j2,2)-c(k2,2))^2 + (a_lab(i2,j2,3)-c(k2,3))^2) + m/s*sqrt((i2-c(k2,4))^2 + (j2-c(k2,5))^2);
                   dis = [dis;d k2];
               end
           end
           if isempty(dis)
           else
           [mind,I] = min(dis(:,1));
           o(i2,j2) = dis(I,2);
           end
       end
   end
   
   for i3 = 1:size(loc,1)
       [row,col] = find(o==i3);
       if isempty(row) && isempty(col)
       else
       rowmean = round(mean(row));
       colmean = round(mean(col));
       c(i3,:)=[a_lab(rowmean,colmean,1) a_lab(rowmean,colmean,2) a_lab(rowmean,colmean,3) rowmean colmean];
       end
   end
   
   iter = iter +1;
    x = iter/(n);
    waitbar(x,f)

   for i4 = 1:size(a,1)
    for j4 = 1:size(a,2)
        for k4 = 1:3
        if o(i4,j4)~=0
        out(i4,j4,k4) = c(o(i4,j4),k4);
        end
    end
    end
end
   
out1 = lab2rgb(out)*255;
figure;
imshow(uint8(out1));

outvid(:,:,:,iter) = uint8(out1);
end

close(f)

%%
for i4 = 1:size(a,1)
    for j4 = 1:size(a,2)
        for k4 = 1:3
        if o(i4,j4)~=0
        out(i4,j4,k4) = c(o(i4,j4),k4);
        end
        end
    end
end
% 
% cform = makecform('lab2srgb');
% out1 = applycform(out,cform);    
out1 = lab2rgb(out)*255;
imshow(uint8(out1));
% %% Edges on the image
d = double(edge((rgb2gray(uint8(out1))),'canny'));
d(find(d==1)) = 255;
d(find(d==0)) = 1;
d(find(d==255)) = 0;
f1 = out1.*d;
figure;
imshow(uint8(f1))
%
[n_rows, n_cols] = size(o);

% Initialize compactness to 0
compactness = 0;

% Loop over each superpixel label in matrix o
for k = 1:max(o(:))
    % Get indices of pixels belonging to superpixel k
    [rows, cols] = find(o == k);
    
    % Calculate centroid of superpixel
    centroid_row = mean(rows);
    centroid_col = mean(cols);
    
    % Calculate average distance of each pixel in the superpixel to centroid
    distances = sqrt((rows - centroid_row).^2 + (cols - centroid_col).^2);
    avg_distance = mean(distances);
    
    % Calculate compactness for superpixel k and add to overall compactness
    s_k = length(rows);
    compactness = compactness + (avg_distance / s_k);
end

disp(compactness);
toc

%%%%%%%%%%%%%%%%%%%%% Part 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3D
load mri;
D = squeeze(D);
A = ind2gray(D,map);

slic = [];
slic0 = [];
meanArr = [];
AllmeanArr = [];

%Loop over different numbers of superpixels with slic method
for i = 10:5:50
    [L,N] = superpixels3(A,i,Method = 'slic');
    imSize = size(A);
    imPlusBoundaries = zeros(imSize(1),imSize(2),3,imSize(3),'uint8');
    for plane = 1:imSize(3)
      BW = boundarymask(L(:, :, plane));
      imPlusBoundaries(:, :, :, plane) = imoverlay(A(:, :, plane), BW, 'green');
    end
    slic=[slic, imPlusBoundaries];
end

%Loop over different numbers of superpixels with slic0 method
for i = 10:5:50
    [L,N] = superpixels3(A,i,Method = 'slic0');
    imSize = size(A);
    imPlusBoundaries = zeros(imSize(1),imSize(2),3,imSize(3),'uint8');
    for plane = 1:imSize(3)
      BW = boundarymask(L(:, :, plane));
      imPlusBoundaries(:, :, :, plane) = imoverlay(A(:, :, plane), BW, 'green');
    end
    slic0=[slic0, imPlusBoundaries];
end

for i = 10:5:50
    meanArr =[];
    for c = [0.0001, 0.001, 0.01, 0.1, 1]
        [L,N] = superpixels3(A,i, method = 'slic0',  Compactness = c);
        pixelIdxList = label2idx(L);
        meanA = zeros(size(A),'like',D);
        for superpixel = 1:N
             memberPixelIdx = pixelIdxList{superpixel};
             meanA(memberPixelIdx) = mean(A(memberPixelIdx));
        end
        meanArr = [meanArr;meanA];
    end
    AllmeanArr = [AllmeanArr,meanArr];
end

implay(A,5)%original mri
implay([slic0;slic],5)% slic0 & slic
implay(AllmeanArr,5)%All results of slic0
