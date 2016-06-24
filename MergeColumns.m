function [new] = MergeColumns(data,mergecols)
[row_m,col_m]=size(mergecols);
row=size(data,1);
k=1; new=[];
for i=1:row
    for r=1:row_m
        k=1;
        for c=1:col_m
            N(k,r)=data(i,mergecols(r,c));
            k=k+1;
        end
    end
    new=[new;N];
end