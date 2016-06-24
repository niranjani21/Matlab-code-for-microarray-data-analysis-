function WriteDataintoFile(data,outfile,outsepa)
a=convertallfiletocell(data);
dlmcell(outfile,a,outsepa);