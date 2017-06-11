function in = compCellArray(cell1,cell2)
% COMPCELLARRAY returns true if each element in the cell array cell1 is
% contained in the cell array cell2
in = all(any(strcmp(repmat(cell1(:),1,numel(cell2)), ...
   repmat(cell2(:)',numel(cell1),1)),2),1);