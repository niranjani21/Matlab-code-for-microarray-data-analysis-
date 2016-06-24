function [aval] = PartofCellString2Double(celldata,cols)
a=celldata;
row_a=size(a,1); 
aval=zeros(row_a,1);
for i=1:row_a
    for j=1:size(cols,2)
        y=ischar(a{i,cols(j)});
        if(y==1),g=str2double(a{i,cols(j)}); aval(i,j)=g; else aval(i,j)=(a{i,cols(j)}); end
    end
end