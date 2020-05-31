I = imread('image.tiff');

BW = imbinarize(I,'adaptive');

figure
imshowpair(I,BW,'montage')

% generate binary image
nx = 1000;
ny = 750;
image_binary = BW;

% specify image domain
x = linspace(-1,1,nx);
y = -linspace(-1,1,ny);

% pad image with zeros in order to enable border at image boundaries
temp = zeros(size(image_binary)+2);
temp(2:end-1,2:end-1) = image_binary;
image_binary = temp;
x = [x(1)-(x(2)-x(1)), x, x(end)+(x(2)-x(1))];
y = [y(1)-(y(2)-y(1)), y, y(end)+(y(2)-y(1))];
[X,Y] = meshgrid(x,y);

% generate edge of the image (subtract eroded image from original image)
image_binary_edge = image_binary-imerode(image_binary,strel('disk',1));

% remove pixels with only one neighbour
image_binary_edge_filtered = imfilter(image_binary_edge,ones(3,3),'same');
image_binary_edge(image_binary_edge_filtered==2) = 0;

% calculate all connected components in image_binary_edge
cc = bwconncomp(image_binary_edge,8);

% initialize vectors for the delaunayTriangulation function 
x_coor = [];
y_coor = [];
constraints = [];
max_dist = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);

% loop over all components

for ii=1:cc.NumObjects
    current = cc.PixelIdxList{ii};
    x_coor_current = X(current);
    y_coor_current = Y(current);

      % reorder coordinates such that they are ordered in a clockwise fashion
      x_coor_reordered = zeros(size(x_coor_current));
      y_coor_reordered = zeros(size(y_coor_current));
      x_coor_reordered(1) = x_coor_current(1);
      y_coor_reordered(1) = y_coor_current(1);
      x_coor_current(1) = [];
      y_coor_current(1) = [];
      kk=2;
      while ~isempty(x_coor_current)
          [index,dist] = knnsearch([x_coor_current,y_coor_current],[x_coor_reordered(kk-1),y_coor_reordered(kk-1)]);

          % if dist is to large, than the current pixel is no neighbouring
          % pixel, this is why we do not at these pixels to the reordered
          % vectors
          if(dist>2*max_dist)
              x_coor_current(index) = [];
              y_coor_current(index) = [];
          else
              x_coor_reordered(kk) = x_coor_current(index);
              y_coor_reordered(kk) = y_coor_current(index);
              x_coor_current(index) = [];
              y_coor_current(index) = [];
              kk = kk + 1;
          end
      end
      x_coor_reordered = x_coor_reordered(1:kk-1); % remove zero entries
      y_coor_reordered = y_coor_reordered(1:kk-1); % remove zero entries

      % take only half of all border samples (this prevents oversampling of 
      % the border)
      x_coor_reordered = x_coor_reordered(1:2:end);
      y_coor_reordered = y_coor_reordered(1:2:end);

      x_coor = [x_coor;x_coor_reordered];
      y_coor = [y_coor;y_coor_reordered];
      constraints_temp = [[length(constraints)+1:length(constraints)+length(x_coor_reordered)]',...
          circshift([length(constraints)+1:length(constraints)+length(x_coor_reordered)]',-1)];
      constraints = [constraints;constraints_temp];
  end

% construct delaunay triangulation
dt = delaunayTriangulation(x_coor,y_coor,constraints);

% maintain only the interior
inside = dt.isInterior();

% Construct a triangulation that represents interior
tr = triangulation(dt(inside, :), dt.Points);

% at the moment, all vertices lie on the edge of the binary image,
% therefore, sample vertices inside the binary image as well:
pointstemp = tr.Points;
connectivityListtemp = tr.ConnectivityList;
pointsinside = zeros(size(X));
for t = 1:size(connectivityListtemp,1)
    vertsXY = pointstemp(connectivityListtemp(t,:),:);
    pointsinside = pointsinside | inpolygon(X,Y, vertsXY(:,1), vertsXY(:,2));
end
pointsinside(1:5:end,:) = 0;
pointsinside(2:5:end,:) = 0;
pointsinside(3:5:end,:) = 0;
pointsinside(4:5:end,:) = 0;
pointsinside(:,1:5:end) = 0;
pointsinside(:,2:5:end) = 0;
pointsinside(:,3:5:end) = 0;
pointsinside(:,4:5:end) = 0;

% construct the triangulation again
dt = delaunayTriangulation([x_coor;X(pointsinside==1)],[y_coor;Y(pointsinside==1)],constraints);
inside = dt.isInterior();
tr = triangulation(dt(inside, :), dt.Points);

% remove points which do not belong to triangle
Points = tr.Points;
ConnectivityList = tr.ConnectivityList;
ii=1;
while(ii<=length(Points))
    if(~isempty(find(ConnectivityList == ii,1)))
        ii = ii + 1;
    else
        Points(ii,:) = [];
        ConnectivityList(ConnectivityList>ii) = ConnectivityList(ConnectivityList>ii)-1;
    end
end
tr = triangulation(ConnectivityList,Points);

% plot the result
figure();
subplot(1,2,1)
imshow(image_binary,[])
title('Binary Image')
subplot(1,2,2)
triplot(tr.ConnectivityList,tr.Points(:,1),tr.Points(:,2))
title('triangulation')
