for ii = 1:3
    
    L{1} = data{ii}(:,61:120)-data{ii}(:,1:60);
    L{2} = data{ii}(:,121:180)-data{ii}(:,1:60);
     for jj = 1:2
         [r p] = corr(gi50',L{jj}');
         [psort pix] = sort(p,'ascend');
         D{ii,jj} = [pws(pix) num2cell([r(pix)' p(pix)' L{jj}(pix,:)])];
     end
end
        