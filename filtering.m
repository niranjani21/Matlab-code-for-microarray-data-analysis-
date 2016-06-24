%function filtering(filename,signalcol_inarray)

a=importdata('normalise_trail.txt'); 
signalcol_inarray=[6 7 8];
% x=iscell(a);
% if x==1
%     a=cell2mat(a);
% end
array=signalcol_inarray;
[r,width_array]=size(array);

for i=1:width_array
    array_mean(i)=mean(a(:,array(i)));
    array_std(i)=std(a(:,array(i)));
end


k=1;
c=[];
for i=1:size(a,1)
    cnt=0;
    for j=1:width_array
        if a(i,array(j)) < array_mean(j)
            b=array_mean(j)-a(i,array(j));
        else
            b=a(i,array(j))-array_mean(j);
        end
        
        if b > array_std(j)
           cnt=cnt+1;
        end
    end

    if cnt == 0
       c(k,:)=a(i,:);
       k=k+1;
    end
end

%end
