function NormalizeAnotateExonArrayData(filename,gc_col,signalcols_inarray,trid_col,outfile)
file=importdata(filename);
data=convertallfiletocell(file);
sig=signalcols_inarray;
[row_sig,col_sig] = size(sig);
[row_a,col_a] = size(data);

norm_data = scale_normalisation(data,sig);
sort_data = sort_myown(norm_data,trid_col,'ascending');
bg_data = ExtractColumnWithFixedElement(sort_data,0,trid_col);
sort_bg = sort_myown(bg_data,gc_col,'ascending');
avg_bg = AverageOverColumnsWithRepeatingEntries(sort_bg,gc_col,sig);


gc=PartofCellString2Double(data,gc_col); 
signal=PartofCellString2Double(data,sig);

for i=1:row_a
    for I=1:size(avg_bg,1)
        if(gc(i,1)==avg_bg(I,1))
            for j=1:col_a
                for Col=1:col_sig
                     if(sig(Col)==j), cnt=1; J=Col; break; else cnt=0; end
                end
                if(cnt==1)
                    if(avg_bg(I,J)>0), data{i,j}=signal(i,J)/avg_bg(I,J+1); end
                end
            end
            break;
        end
    end
end
dlmcell(outfile,data,',');