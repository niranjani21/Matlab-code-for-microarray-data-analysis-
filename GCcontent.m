function GCcontent(file,GC_col,outfile,outsepa)
b=importdata(file);
a=convertallfiletocell(b);
row_a=size(a,1);

for i=1:row_a
    base = basecount(a{i,GC_col});
    B=struct2cell(base);
    GC=(B{2,1}+B{3,1})*100/(B{1,1}+B{2,1}+B{3,1}+B{4,1});
    a{i,GC_col}=GC;  
end
dlmcell(outfile,a,outsepa);