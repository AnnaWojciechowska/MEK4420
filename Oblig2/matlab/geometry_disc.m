% box of draught D and length (width) L, with D as unit length
% left vertical side
depth = 1; %Depth
length = 2; % Length

n_side_segments = 10; % number of segments on side
n_bottom_segments = 20; % number of segments on botttom
n_total_segmnets = 2* n_side_segments + n_bottom_segments;

dy = depth/n_side_segments;
dx = length/n_bottom_segments;

% each is segment is  marked as :M-P
% (xm[i],ym[i]) ------------------------- (xp[i],yp[i]) 

%prealocating arrays for performance reasons:
x_m = zeros(1, n_total_segmnets);
y_m = zeros(1, n_total_segmnets);
x_p = zeros(1, n_total_segmnets);
y_p = zeros(1,n_total_segmnets);


% initializing coordinates:
x_m(1) = -length/2; y_m(1) = 0; %since -dy*(1-1) is 0
x_p(1) = -length/2; y_p(1) = -dy;
% coord is 2d array of segment coordinates
% shape is 4 x n_total_segmnets
% each row is one segment with its M and P point coord
coord=[x_m(1),y_m(1),x_p(1),y_p(1)];

% first we initialize coordinates array by 1st side
for i=2 : n_side_segments
    x_m(i) = -length/2;  y_m(i) = -dy*(i-1);
    x_p(i) = -length/2;  y_p(i) = -dy*i;
    segment=[x_m(i),y_m(i),x_p(i),y_p(i)];
    coord=[coord;segment];

end


% now we caluclate points coordinates allong bottom
for i= 1+n_side_segments : n_side_segments+n_bottom_segments
    i1=i-n_side_segments;
    x_m(i)=-length/2+dx*(i1-1); y_m(i)=-depth;
    x_p(i)=-length/2+dx*i1;  y_p(i)=-depth;
    segment=[x_m(i),y_m(i),x_p(i),y_p(i)];
    coord=[coord;segment];
   
end



% and the last side
for i=1+n_side_segments+n_bottom_segments:n_side_segments+n_bottom_segments+n_side_segments;
    i1=i-n_side_segments-n_bottom_segments;
    x_m(i)=length/2; y_m(i)=-depth+dy*(i1-1);
    x_p(i)=length/2; y_p(i)=-depth+dy*i1;
    segment=[x_m(i),y_m(i),x_p(i),y_p(i)];
    coord=[coord;segment];
end



% testing part to compare if same coords are generated in python script
fileID = fopen('matlab_points.txt', 'w');
fprintf(fileID, '%d ', x_m);
fprintf(fileID, '\n');
fprintf(fileID, '%d ', y_m);
fprintf(fileID, '\n');
fprintf(fileID, '%d ', x_p);
fprintf(fileID, '\n');
fprintf(fileID, '%d ', y_p);
fclose(fileID);


% saving coord to box.dat that is later used in another script
save -ascii box.dat coord;

%clean the plot
clf

% plot geometry
% since segments are x_p(i)... x_m(i) to get all points the last
% x_p(n_total_segmnets) must be added at end of array
x_plot = [x_m,x_p(n_total_segmnets)];
y_plot = [y_m,y_p(n_total_segmnets)];

hold on

axis ([-1.01 1.01 -1.1 0.01])
handle = plot(x_plot,y_plot, 'k *');
%axis off
get(handle);
set(handle,'MarkerSize',[10])
set(gca,'FontSize',20)