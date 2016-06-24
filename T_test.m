function T_test(filename,colsToBeTested_inarray)
a=importdata(filename);
[ra,widthA]=size(a);
array=colsToBeTested_inarray;
[warray,width_array]=size(array);
for i=1:length(array)
    x=array(i,1);
    y=array(i,2);
    [h,p] = mattest(a(:,x),a(:,y))
end

 