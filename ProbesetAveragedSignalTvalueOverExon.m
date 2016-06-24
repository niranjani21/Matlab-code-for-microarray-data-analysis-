function ProbesetAveragedSignalTvalueOverExon(file,acols,sgnlcols,mergecols,ccols,Pval,fccutoff,log_base,shift,outfile)
a=importdata(file);
data=convertallfiletocell(a);
gene_sig=AveragedGeneSignalsOverTranscript(data,acols,sgnlcols,mergecols,log_base,shift);
sort_data=sort_myown(data,acols(2),'ascending');
disp(['data sorted']);
signal = PartofCellString2Double(sort_data,sgnlcols);
disp(['signal columns retrieved']);
ids = PartofCellString2Double(sort_data,acols);
[row,col]=size(data);
[row_m,col_m]=size(mergecols);
[row_c,col_c]= size(ccols);
[row_g,col_g]=size(gene_sig);
for i=1:row_m
    lcols(1,i)=i;
end

exon_count=0;
last_exon=0; k=1; tv=[];
for i=1:row
    tv=[tv;signal(i,:)];
    if(i<row)
        if(ids(i,2) ~= ids(i+1,2))
            last_exon=0;
            if(ids(i,1) ~= ids(i+1,1))
                last_exon=1;
                exon_count=exon_count+1;
            else
                exon_count=exon_count+1;
            end
            merge_data=MergeColumns(tv,mergecols);
            for L=1:row_c
                fc(L,1)=FoldchangeCalculation(merge_data,ccols(L,1),ccols(L,2));
            end
            lv=log_transformOfDoubleData(merge_data,lcols,log_base,shift);
            for L=1:row_c
                tval(L,1)=TvalueCalculation(lv,ccols(L,1),ccols(L,2));
            end
            V=size(data,1)-1;
            tab_tval = tinv(Pval,V);
            cnt=0;
            for L=1:row_c
                if(tab_tval < tval(L,1) && fc(L,1) > fccutoff)
                    cnt=cnt+1;
                end
            end

            if(cnt > 1)
               sda(k,1)=ids(i,1);
               sda(k,2)=ids(i,2);
               for j=1:row_c
                   sda(k,2+j)=fc(j,1);
                   sda(k,2+row_c+j)=tval(j,1);

               end
               for G=1:col_g
                   if(gene_sig(G,1) == sda(k,1))
                       for C=1:row_c 
                           sda(k,2+2*row_c+C) = (mean(lv(:,ccols(C,1)))/gene_sig(G,3+ccols(C,1))) - (mean(lv(:,ccols(C,2)))/gene_sig(G,3+ccols(C,2)));    
                       end
                       break;
                   end
               end
               
               k=k+1;
            end
            tv=[];
        end
    else
        merge_data=MergeColumns(tv,mergecols);
        for L=1:row_c
            fc(L,1)=FoldchangeCalculation(merge_data,ccols(L,1),ccols(L,2));
        end
        lv=log_transformOfDoubleData(merge_data,lcols,log_base,shift);
        for L=1:row_c
            tval(L,1)=TvalueCalculation(lv,ccols(L,1),ccols(L,2));
        end
        V=size(data,1)-1;
        tab_tval = tinv(Pval,V);
        cnt=0;
        for L=1:row_c
            if(tab_tval < tval(L,1) && fc(L,1) > fccutoff)
                cnt=cnt+1;
            end
        end

        if(cnt > 1)
           sda(k,1)=ids(i,1);
           sda(k,2)=ids(i,2);
           for j=1:row_c
               sda(k,2+j)=fc(j,1);
               sda(k,2+row_c+j)=tval(j,1);            
           end
           for G=1:col_g
               if(gene_sig(G,1) == sda(k,1))
                   for C=1:row_c                    
                       sda(k,2+2*row_c+C) = (mean(lv(:,ccols(C,1)))/gene_sig(G,3+ccols(C,1))) - (mean(lv(:,ccols(C,2)))/gene_sig(G,3+ccols(C,2)));
                   end
                   break;
               end
           end
        end
    end
    
    if(last_exon==1), exon_count=0; end
end    
%dlmwrite(outfile,sda,',');