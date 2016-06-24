function [out] = AveragedGeneSignalsOverTranscript(data,acols,sgnlcols,mergecols,log_base,shift)
%a=importdata('out.txt'); 
%data=convertallfiletocell(a); acols=[2 3]; sgnlcols=[7 8 9 10 11 12 13 14 15]; mergecols=[1 2 3;4 5 6;7 8 9]; log_base=2; shift=1;
sort_data=sort_myown(data,acols(2),'ascending');
signal=PartofCellString2Double(sort_data,sgnlcols);
ids=PartofCellString2Double(sort_data,acols);
col_sig=size(sgnlcols,2);
[row,col] = size(data);
exon_count=0; tv=[]; k=1; 
for i=1:size(mergecols,1)
    lcols(1,i)=i;
end
sda=zeros(1,3+size(lcols,2));
for i=1:row
    if(i<row)
        tv=[tv;signal(i,:)];
 
        if(ids(i,2) ~= ids(i+1,2))
            tvt=MergeColumns(tv,mergecols);
            sda(k,3)=sda(k,3)+size(tvt,1);
            for j=1:size(lcols,2)
                 lv=log_transformOfDoubleData(tvt(:,j),lcols,log_base,shift);
                 sda(k,3+j)=sda(k,3+j)+mean(lv);
            end
            tv=[]; exon_count=exon_count+1;

            if(ids(i,1) ~= ids(i+1,1))             
                sda(k,1)=ids(i,1);
                sda(k,2)=exon_count;
                for j=1:size(lcols,2)
                    sda(k,3+j)=sda(k,3+j)/exon_count;
                end
                k=k+1;
                exon_count=0;
                sda(k,:)=0;
            end

        end
    else
        if(size(tv,1) ~= 0 )
            tvt=MergeColumns(tv,mergecols);
            sda(k,3)=sda(k,3)+size(tvt,1);
            for j=1:size(lcols,2)
                lv=log_transformOfDoubleData(tvt(:,j),lcols,log_base,shift);
                sda(k,3+j)=(sda(k,3+j)+mean(lv/exon_count));
            end
            sda(k,1)=ids(i,1);
            sda(k,2)=exon_count;
        end        
    end
    
end
[R,C]=size(sda); cnt=0;
for j=1:C
    if(sda(R,j) == 0), cnt=cnt+1; end
end
if(cnt ~= C),out=sda; else out=sda(1:R-1,:); end