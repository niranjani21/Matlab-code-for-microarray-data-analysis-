function [out] = AverageOverColumnsWithRepeatingEntries(data,ref_col,signl_col)
%d=importdata('out.txt'); ref_col=6; signl_col=[7 8];
data=convertallfiletocell(data);
[a] = sorting(data,ref_col,'ascending');
[row_a,col_a] = size(a);
col_sig=size(signl_col,2);

for i=1:row_a
    k=1;   
    if(ischar(a{i,ref_col})==1),ref(i,1)=str2double(a{i,ref_col}); else ref(i,1)=a{i,ref_col}; end
    for j=signl_col(1):signl_col(col_sig)
        y=ischar(a{i,j});
        if(y==1),aval(i,k)=str2double(a{i,j}); else aval(i,k)=a{i,j}; end
        k=k+1;
    end
end

tv=[]; count=0; i=1; k=1;
while(i <= row_a)
   
    for w=i:row_a
        tv=[tv;aval(i,:)];
        count=count+1; 
        if(i<row_a)
            if(ref(i,1) ~= ref(i+1,1))
                sda(k,1)=ref(i,1);
                for j=1:col_sig
                    sda(k,j+1)=sum(tv(:,j))/count;
                end
                tv=[]; count=0; k=k+1; i=i+1; break;
            end
        else
            sda(k,1)=ref(i,1);
            for j=1:col_sig
                sda(k,j+1)=sum(tv(:,j))/count;
            end
        end
        i=i+1;
    end
end
out=sda;
