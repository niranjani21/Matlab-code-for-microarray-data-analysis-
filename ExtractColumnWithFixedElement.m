function [out] = ExtractColumnWithFixedElement(data,fixed_element,e_col)
%data=sort_data; fixed_element=0; e_col=2;
a=convertallfiletocell(data);
[row_a,col_a]=size(a);
k=1;
for i=1:row_a
    aval=PartofCellString2Double(a,e_col);
    if(aval(i,1)==fixed_element)
        for j=1:col_a
            b{k,j}=a{i,j};  
        end
        k=k+1;
    end
end
out=b;