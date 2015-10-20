%hazel's change
data=[1 3 3;3 4 5];

f = figure;
colnames={'vol1','vol2','vol3'};
rownames={'jij','ee'};
t=uitable(f,'Data',data,'ColumnName',colnames,'Rowname',rownames);
t.Position(1) = t.Extent(1);
t.Position(2) = t.Extent(2);

     