function [c] = log_transform(data,cols_inarray,log_base,shift)

[row_a,col_a]=size(data);
sig=cols_inarray;
a=data;
for i=1:row_a
    for j=1:col_a
        cnt=0;
        if(ismember(j,sig)==1),cnt=1; end
        if(cnt == 1)
            y=ischar(a{i,j});
            if(y==1),aval=str2double(a{i,j}); else aval=a{i,j}; end
            c{i,j}=log2(aval+shift)/log2(log_base);
        else
            c{i,j}=a{i,j};
        end
            
    end
end
