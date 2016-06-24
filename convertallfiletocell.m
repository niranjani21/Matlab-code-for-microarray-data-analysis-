function [cell_file]=convertallfiletocell(data)
%data=importdata('GC_trail.txt');
a=data;

x=isstruct(a);
if x==1
    [textA_row,textA_col]=size(a.textdata);
    [dataA_row,dataA_col]=size(a.data);
    for i=1:dataA_row
        for j=1:dataA_col
            d{i,j}=num2str(a.data(i,j));
        end
    end
    
    A=a.textdata;
    k=textA_col+1;
   

for j=1:dataA_col
    for i=1:size(d,1)
        A{i,k}=d{i,j};
    end
    k=k+1;
end
else
    A=num2cell(a);
end

y=iscell(a);
if y==1
    A=a;
end
cell_file=A;


