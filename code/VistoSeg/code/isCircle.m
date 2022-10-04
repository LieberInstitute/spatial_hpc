
nSpots = size(tbl, 1);
disp([num2str(nSpots),' Visium spots detected'])
    

for i = 1:nSpots 

    crow = round(table2array(tbl(:, 5)));
    ccol = round(table2array(tbl(:, 6)));
    idx = mask(crow(i), ccol(i));
    spot = regionprops(mask==idx);
    
    for C = 1:numel(O)
        signal = struct2table(regionprops(mask==idx & BW.(O{C})>0));
        points = signal.Centroid;        
        isincircle = sum((points - [ccol(i) crow(i)]).^2,2)<= R^2;

        %check
%         [tempx,tempy] = find(mask == idx);
%         temp = BW.(O{C})(min(tempx):max(tempx),min(tempy):max(tempy));
%         imshow(temp)
%         viscircles(size(temp)/2, repelem(R, 1), 'Color', 'r');

        count.(O{C})(i) = length(find(isincircle));
        prop.(O{C})(i) = sum([signal.Area])/spot.Area;
    end
      
      if mod(i,100) == 0
        disp([num2str(i),' spots finished in time ', num2str(toc),'s'])
      end

end
