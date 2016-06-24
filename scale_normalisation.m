function [c] = scale_normalisation(data,signalcols_inarray)
%data=importdata('GC_trail2.txt'); signalcols_inarray=[7 8 9];
a=convertallfiletocell(data);

[row_a,col_a]=size(a);
sig=signalcols_inarray;
[row_sig,col_sig]=size(sig);

k=1;
for j=sig(1):sig(col_sig)
    for i=1:row_a
        y=ischar(a{i,j});
        if(y==1), g=str2double(a{i,j}); aval(i,k)=g; else aval(i,k)=a{i,j}; end
    end
    sum_total(k)=sum(aval(:,k));
    k=k+1;
end

sum_mean=mean(sum_total);

for j=1:col_a
    cnt=0;
    for J=1:col_sig
        if(j==sig(J)), cnt =1; s_tot=sum_total(J); break; end
    end                
    if(cnt==1)
        for i=1:row_a , c{i,j}=aval(i,J)*(sum_mean/s_tot); end
    else
        for i=1:row_a , c{i,j}=a{i,j}; end
    end

end


