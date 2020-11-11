img = imread('1.bmp');

%znajdz źrenicę
bin = im2bw(img,0.1);
[centers,radii] = imfindcircles(bin, [30,1000], 'ObjectPolarity','dark')

%wytnij z obrazu pierścień o wymiarach (r_źrenicy; r_źrenicy*3)
inner_rad = radii;
outer_rad = radii*3;

center_x = centers(1);
center_y = centers(2);

img_size = size(img);

[x,y] = meshgrid(1:img_size(2),1:img_size(1));

distance = (x-center_x).^2+(y-center_y).^2;
mask = uint8(distance<outer_rad^2 & distance>inner_rad^2);

mask(mask==1)=255;

masked_img = bitand(img,mask);

resized_img = masked_img(uint16(centers(2) - outer_rad):uint16(centers(2) + outer_rad),uint16(centers(1)-outer_rad):uint16(centers(1) + outer_rad));
imshow(resized_img);

% wyszukaj powieki
filtered_img = edge(resized_img,'canny');

imshow(filtered_img);
img_size = size(filtered_img)

[H,theta,rho] = hough(filtered_img);
P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(filtered_img,theta,rho,P,'FillGap',img_size(1)/2,'MinLength',7);
figure, imshow(filtered_img), hold on

clostest_top = [[0,0]; [0,0]];
clostest_bot = [img_size; img_size];
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   % bottom
   if xy(1,2) > img_size(1)/2 && xy(2,2) > img_size(1)/2
       if (xy(1,2) + xy(2,2)) / 2 < (clostest_bot(1,2) + clostest_bot(2,2)) / 2
           clostest_bot = [[xy(1,1), xy(1,2)]; [xy(2,1), xy(2,2)]];
       end
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
   end
   
   % top
   if xy(1,2) < img_size(1)/2 && xy(2,2) < img_size(1)/2
       if (xy(1,2) + xy(2,2)) / 2 > (clostest_top(1,2) + clostest_top(2,2)) / 2
           clostest_top = [[xy(1,1), xy(1,2)]; [xy(2,1), xy(2,2)]];
       end
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
   end
end

% usuń powieki
if clostest_bot ~= [img_size; img_size]
    plot(clostest_bot(:,1),clostest_bot(:,2),'LineWidth',2,'Color','red');
    
    coefficients = polyfit([clostest_bot(1,1), clostest_bot(2,1)], [clostest_bot(1,2), clostest_bot(2,2)], 1);
    a = coefficients (1);
    b = coefficients (2);
    
    [x,y] = meshgrid(1:img_size(2),1:img_size(1));
    mask = uint8(a * x + b > y);
    mask(mask==1)=255;
    resized_img = bitand(resized_img,mask);
end

if clostest_top ~= [img_size; img_size]
    plot(clostest_top(:,1),clostest_top(:,2),'LineWidth',2,'Color','red');
    
    coefficients = polyfit([clostest_top(1,1), clostest_top(2,1)], [clostest_top(1,2), clostest_top(2,2)], 1);
    a = coefficients (1);
    b = coefficients (2);
    
    [x,y] = meshgrid(1:img_size(2),1:img_size(1));
    mask = uint8(a * x + b < y);
    mask(mask==1)=255;
    resized_img = bitand(resized_img,mask);
end
figure, imshow(resized_img)

% wyszukaj tęczówkę
binared = imbinarize(resized_img,"adaptive");
figure, imshow(binared)

edge_left = 5
edge_right = img_size(1) - 5
for i = 5:img_size(2)
    if binared(uint8(img_size(1)/2),i) == 0
        edge_left = i
        break
    end
end

for i = img_size(2)-5:1
    if binared(uint8(img_size(1)/2),i) == 0
        edge_right = i
        break
    end
end

% zostaw tylko tęczówkę
outer_r = min(img_size(2)/2 - edge_left,edge_right - img_size(2)/2)
outer_r = outer_r-15
center_x = img_size(1)/2
center_y = img_size(2)/2

[x,y] = meshgrid(1:img_size(2),1:img_size(1));

distance = (x-center_x).^2+(y-center_y).^2;
mask = uint8(distance<outer_r^2 & distance>inner_rad^2);

mask(mask==1)=255;

resized_img = bitand(resized_img,mask);
figure;imshow(resized_img)

% przekształć na współrzędne biegunowe
L=floor(2*pi*outer_r);

t=[45:360/L:404+1-360/L];   
R_inner=inner_rad;x_inner=R_inner*sind(t)+center_x;y_inner=R_inner*cosd(t)+center_y; 
hL2=plot(x_inner,y_inner,'m');
x_ref_inner=hL2.XData;y_ref_inner=hL2.YData;
                                 
R=outer_r;x=R*sind(t)+center_x;y=R*cosd(t)+center_y; 
hL1=plot(x,y,'m'); % axis equal;grid on;  

x_ref=hL1.XData;y_ref=hL1.YData;

Sx={};Sy={};
for k=1:1:numel(hL1.XData)
    x_ref(k)
    center_x+x_ref_inner(k)
    Lx=floor(linspace(x_ref(k),x_ref_inner(k),ceil(R-inner_rad)));
    Ly=floor(linspace(y_ref(k),y_ref_inner(k),ceil(R-inner_rad)));
    Sx=[Sx Lx'];Sy=[Sy Ly'];
end

sx=cell2mat(Sx);sy=cell2mat(Sy);
[s1 s2]=size(sx);

B1=uint8(zeros(s1,s2));

for n=1:1:s2
    for k=1:1:s1
        B1(k,n)=resized_img(sx(k,n),sy(k,n)); 
    end
end
img_teczowka=uint8(zeros(s1,s2,1));
img_teczowka(:,:,1)=B1;
figure;imshow(img_teczowka);