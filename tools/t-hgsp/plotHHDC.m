function [h] = plotHHDC(hhdcIn, valor)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
alpha = 0.1;
% hhdcIn(hhdcIn == 0) = NaN;
hhdcIn = padarray(padarray(hhdcIn,1,'pre'),1,'post');
hhdcIn = permute(hhdcIn,[2 1 3]);
hhdcIn = padarray(padarray(hhdcIn,1,'pre'),1,'post');
hhdcIn = permute(hhdcIn,[2 1 3]);
hhdcIn = permute(hhdcIn,[3 2 1]);
hhdcIn = padarray(padarray(hhdcIn,1,'pre'),1,'post');
hhdcIn = permute(hhdcIn,[3 2 1]);
subcubo = hhdcIn;
if (valor == 1)
    subcubo = flip(hhdcIn,3);
end
aux = max(subcubo,[],3);
% subcubo = subcubo ./ aux;
n = 1;
for i = 1:size(subcubo,3)
    TheImage = flip(subcubo(:,:,size(subcubo,3)+1-i));
    [M, N, P] = size(TheImage);
    [Y,X] = ndgrid(1:M,1:N);
    Z = zeros(size(X)) + n;
    amap = ones(size(TheImage));
    if(i == 2 || i == size(subcubo,3))
        amap(TheImage == 0.0) = alpha; % This value changes the alpha transparency
    else
        amap(TheImage == 0.0) = 0; % This value changes the alpha transparency
    end
    amap(:,1) = 0;
    amap(:,end) = 0;
    amap(1,:) = 0;
    amap(end,:) = 0;
    
    surf(X, Y, Z, TheImage, 'edgecolor','none','AlphaData',amap,'FaceAlpha','flat','AlphaDataMapping','none');
    % I usually change the X and Y lims so it looks better
    set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
    colormap('jet')
    hold on
    n = n + 1;
end
% n = 0;
for i = size(subcubo,3):-1:1
    TheImage = flip(subcubo(:,:,size(subcubo,3)+1-i));
    [M, N, P] = size(TheImage);
    [Y,X] = ndgrid(1:M,1:N);
    Z = zeros(size(X)) + n;
    amap = ones(size(TheImage));
%     if(i == 2)
%         amap(TheImage == 0.0) = alpha; % This value changes the alpha transparency
%     else
        amap(TheImage == 0.0) = 0; % This value changes the alpha transparency
%     end
    surf(X, Y, Z, TheImage, 'edgecolor','none','AlphaData',amap,'FaceAlpha','flat','AlphaDataMapping','none');
    % I usually change the X and Y lims so it looks better
    set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
    colormap('jet')
    hold on
    n = n - 1;
end

view(-135,45) % Angle of view of the cube
daspect([1 1 1]) % Change the 3 value to "compress" the cube
% hhdc =hhdcIn;
% hhdc = 0;
h=0;
end