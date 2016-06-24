function [added_data] = add_lowestsignal(filename,signalcol_inarray)

A=importdata('normalise_trail.txt'); signalcol_inarray=[6 7 8];
a=convertallfiletocell(A);
[row_a,col_a]=size(a);
sig=signalcol_inarray;
[row_sig,col_sig]=size(sig);

for i=1:row_a
    for j=sig(1):sig(col_sig)
        aval(i,j)=a{i,j};
    end
end

M=min(aval(:,sig(1):sig(col_sig)));
MIN=min(M);

add_val=0.02*MIN;

acol=[];
af=1;
j=1;
while af <= col_a
    cnt=0;
    for i=1:col_sig
        if af == sig(i)
            cnt=cnt+1;
        end
    end
    if cnt == 0
        acol(j)=af;
        j=j+1;
    end
    af=af+1;
end

for i=1:row_a
    for j=1:size(acol,1)
        c{i,j}=a{i,acol(j)};
    end
    k=length(acol)+1;
    for j=1:col_sig
        c{i,k}=a{i,sig(j)}+add_val;
        k=k+1;
    end
end
added_data=c;

