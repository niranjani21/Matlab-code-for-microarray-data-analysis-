classdef ma_analysis
   
    methods (Static,  Access = public)
        
        function [added_data] = add_lowestsignal(filename,signalcol_inarray)
            if(exist(filename) == 0)
                disp(filename,'does not exist');
            else
                A=importdata(filename);
                a=convertallfiletocell(A);
                [row_a,col_a]=size(a);
                sig=signalcol_inarray;
                [~,col_sig]=size(sig);
               
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
            end
        end
        
        %********************************************************************************
                
        function [out] = AverageOverColumnsWithRepeatingEntries(data,ref_col,signl_col)
            if(exist(filename) == 0)
                disp(filename,'does not exist');
            else
                data=convertallfiletocell(data);
                [a] = sorting(data,ref_col,'ascending');
                [row_a,~] = size(a);
                col_sig=size(signl_col,2);
                
                for i=1:row_a
                    k=1;   
                    if(ischar(a{i,ref_col})==1),ref(i,1)=str2double(a{i,ref_col}); else ref(i,1)=a{i,ref_col}; end
                    for j=signl_col(1):signl_col(col_sig)
                        y=ischar(a{i,j});
                        if(y==1),aval(i,k)=str2double(a{i,j}); else aval(i,k)=a{i,j}; end
                        k=k+1;
                    end
                end

                tv=[]; count=0; i=1; k=1;
                while(i <= row_a)

                    for w=i:row_a
                        tv=[tv;aval(i,:)];
                        count=count+1; 
                        if(i<row_a)
                            if(ref(i,1) ~= ref(i+1,1))
                                sda(k,1)=ref(i,1);
                                for j=1:col_sig
                                    sda(k,j+1)=sum(tv(:,j))/count;
                                end
                                tv=[]; count=0; k=k+1; i=i+1; break;
                            end
                        else
                            sda(k,1)=ref(i,1);
                            for j=1:col_sig
                                sda(k,j+1)=sum(tv(:,j))/count;
                            end
                        end
                        i=i+1;
                    end
                end
                out=sda;
            end
        end
        
        %********************************************************************************
       
        function [cnt] = CheckForColumn(col_no,check_cols)
            col=size(check_cols,2);
            for i=1:col
                if(check_cols(i)==col_no), cnt(1)=1; cnt(2)=i; break; else cnt(1)=0; cnt(2)=0; end
            end
        end
        
        %********************************************************************************
       
        function [cnt,sc] = CheckForColumnNumber(col_no,check_cols)
            col=size(check_cols,2);
            for i=1:col
                if(check_cols(i)==col_no),sc=i; cnt=1; break; else cnt=0;end
            end
        end
        
        %********************************************************************************
       
        function comparing(filename_a,filename_b,colna_inarray,colnb_inarray,outfile)
            colna=colna_inarray;
            colnb=colnb_inarray;
            [col_colna]=size(colna,2);
            [col_colnb]=size(colnb,2);
            if(exist(filename_a) == 0 || exist(filename_b) == 0 )
                disp('One or many input files missing');
            else
                a=importdata(filename_a);
                b=importdata(filename_b);

                [a]=sort_myown(a,colna(1,1),'ascending');
                [b]=sort_myown(b,colnb(1,1),'ascending');

                [row_a,col_a]=size(a);
                [row_b,col_b]=size(b);

                if(col_colna == col_colnb) 
                    check_count = 0;
                    for k=1:col_colna
                        if (colna(k) > 0 && colna(k) <= col_a )
                            if (colnb(k) >0 && colnb(k) <= col_b)
                                check_count=check_count+1;
                            end
                        end
                    end
                end

                % convert cell to double for comparison
                aval=cell2mat(a);
                
                common_items=0;
                if(check_count == 0)
                    disp(['warning: cannot compare, column numbers are out of range']);
                else
                    i=1; ic=1; b_start=1; 
                     while i <= row_a
                        for j=b_start:row_b
                            count=0;
                            if (aval(i,colna(1)) > bval(j,colnb(1))) 
                                b_start=j;
                            else if aval(i,colna(1)) < bval(j,colnb(1))
                                    break;
                                else     
                                    for cc=1:col_colna
                                        if aval(i,colna(cc)) == bval(j,colnb(cc))
                                            count=count+1;   
                                        end           
                                    end

                                    if count == col_colna
                                       common_items=common_items+1;

                                       for cc=1:col_colna    
                                           c{ic,cc}=a{i,colna(cc)};
                                       end
                                       
                                       al = col_colna+1;
                                       af=1;

                                       while (af <= col_a)
                                           cola_cnt=0;
                                           for cnt=1:col_colna
                                               if (af == colna(cnt))
                                                  cola_cnt=cola_cnt+1;
                                               end
                                           end

                                           if cola_cnt==0
                                              c{ic,al}=a{i,af};
                                              al=al+1;
                                           end
                                           af=af+1;
                                       end

                                       bl=col_a+1;
                                       bf=1;

                                       while (bf <= col_b)
                                           colb_cnt=0;
                                           for cnt=1:col_colnb
                                               if (bf == colnb(cnt))
                                                   colb_cnt=colb_cnt+1;
                                               end
                                           end

                                           if colb_cnt==0
                                              c{ic,bl}=b{j,bf};
                                              bl=bl+1;
                                           end
                                           bf=bf+1;
                                       end
                                        ic=ic+1;
                                    end
                                end
                            end
                        end
                        i=i+1;
                    end
                end

                if(size(c,1)==0)
                    disp('no common items found');
                else
                    disp([num2str(size(c,1)),' common items found, writing into file']);
                    dlmcell(outfile,c,',');
                end
            end
        end
    
        %********************************************************************************
            
        function [cell_file]=convertallfiletocell(filename)
            if(exist(filename) == 0)
                disp(filename,'does not exist');
            else
                a=importdata(filename);
                x=isstruct(a);
                if x==1
                    [~,textA_col]=size(a.textdata);
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
            end
        end
    
        %********************************************************************************
            
        function [out] = ExtractColumnWithFixedElement(data,fixed_element,e_col)
            a=convertallfiletocell(data);
            [row_a,col_a]=size(a);
            k=1;
            for i=1:row_a
                aval=PartofCellString2Double(a,e_col);
                if(aval(i,1)==fixed_element)
                    for j=1:col_a
                        b{k,j}=a{i,j};  
                    end
                    k=k+1;
                end
            end
            out=b;
        end

        %********************************************************************************
        
        function filtering(filename,signalcol_inarray)
            if(exist(filename) == 0)
                disp(filename,'does not exist');
            else
                a=importdata(filename);
            
                array=signalcol_inarray;
                [~,width_array]=size(array);

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
            end
        end
        
        %********************************************************************************
        
        function [fc] = FoldchangeCalculation(filename,col_1,col_2)
            if(exist(filename) == 0)
                disp(filename,'does not exist');
            else
                a=importdata(filename); 
                data=convertallfiletocell(a);

                o_col=data(:,col_1);
                t_col=data(:,col_2);

                o_mean=mean(o_col);
                t_mean=mean(t_col);
                mean_max=max(o_mean,t_mean); mean_min=min(o_mean,t_mean);

                if mean_min > 0
                    if(o_mean > t_mean)
                        fc = mean_max/mean_min;
                    else
                        fc = -1 * (mean_max/mean_min);
                    end
                else
                    fc = 0;
                end
            end
        end

        %********************************************************************************
        
        function GCcontent(filename,GC_col,outfile,outsepa)
            if(exist(filename) == 0)
                disp(filename,'does not exist');
            else
                b=importdata(filename);
                a=convertallfiletocell(b);
                row_a=size(a,1);

                for i=1:row_a
                    base = basecount(a{i,GC_col});
                    B=struct2cell(base);
                    GC=(B{2,1}+B{3,1})*100/(B{1,1}+B{2,1}+B{3,1}+B{4,1});
                    a{i,GC_col}=GC;  
                end
                dlmcell(outfile,a,outsepa);
            end
        end
        
        %********************************************************************************
        
        function [c] = log_transform(data,cols_inarray,log_base,shift)
            [row_a,col_a]=size(data);
            sig=cols_inarray;
            a=data;
            for i=1:row_a
                for j=1:col_a
                    cnt=0;
                    if(ismember(j,sig)==1),cnt=1; end
                    if(cnt == 1)
                        y=ischar(a{i,j});
                        if(y==1),aval=str2double(a{i,j}); else aval=a{i,j}; end
                        c{i,j}=log2(aval+shift)/log2(log_base);
                    else
                        c{i,j}=a{i,j};
                    end

                end
            end
        end

        %********************************************************************************
        
        function [c] = log_transformOfDoubleData(data,cols_inarray,log_base,shift)
            [row_a,col_a]=size(data);
            sig=cols_inarray;
            a=data;
            for i=1:row_a
                for j=1:col_a
                    cnt=0;
                    if(ismember(j,sig)==1),cnt=1; end
                    if(cnt == 1)
                        c(i,j)=log2(a(i,j)+shift)/log2(log_base);
                    else
                        c(i,j)=a(i,j);
                    end

                end
            end
        end
        
        %********************************************************************************
                
        function [new] = MergeColumns(data,mergecols)
            [row_m,col_m]=size(mergecols);
            row=size(data,1);
            k=1; new=[];
            for i=1:row
                for r=1:row_m
                    k=1;
                    for c=1:col_m
                        N(k,r)=data(i,mergecols(r,c));
                        k=k+1;
                    end
                end
                new=[new;N];
            end
        end
        
        %********************************************************************************
        
        function NormalizeAnotateExonArrayData(filename,gc_col,signalcols_inarray,trid_col,outfile)
            if(exist(finlename) == 0)
                disp(filename,'does not exist');
            else
                file=importdata(filename);
                data=convertallfiletocell(file);
                sig=signalcols_inarray;
                [~,col_sig] = size(sig);
                [row_a,col_a] = size(data);

                norm_data = scale_normalisation(data,sig);
                sort_data = sorting(norm_data,trid_col,'ascending');
                bg_data = ExtractColumnWithFixedElement(sort_data,0,trid_col);
                sort_bg = sorting(bg_data,gc_col,'ascending');
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
            end
        end
        
        %********************************************************************************
        
        function [aval] = PartofCellString2Double(celldata,cols)
            a=celldata;
            row_a=size(a,1); 
            aval=zeros(row_a,1);
            for i=1:row_a
                for j=1:size(cols,2)
                    y=ischar(a{i,cols(j)});
                    if(y==1),g=str2double(a{i,cols(j)}); aval(i,j)=g; else aval(i,j)=(a{i,cols(j)}); end
                end
            end
        end
        
        %********************************************************************************
                
        function [c] = scale_normalisation(filename,signalcols_inarray)
            if(exist(finlename) == 0)
                disp(filename,'does not exist');
            else
                data=importdata(filename);
                a=convertallfiletocell(data);

                [row_a,col_a]=size(a);
                sig=signalcols_inarray;
                [~,col_sig]=size(sig);

                k=1;
                for j=sig(1):sig(col_sig)
                    aval(:,k)=cell2mat(a{:,j});
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
            end
        end
        
        %********************************************************************************
         
        function [c]=sorting(data,sort_colno,sorttype)
            [a]=convertallfiletocell(data);

            if(strcmp('ascending',sorttype)==1)
                s=sort_colno;
            else
                S=sprintf('-%d',sort_colno);
                s=str2double(S);
            end
            c=sortrows(a,s);
        end
              
        %********************************************************************************
         
        function [h_vals,p_vals,t_vals] = TvalueCalculation(data,comparecols_inarray)
            [c_rows,c_cols] = size(comparecols_inarray);
            h_vals=zeros(c_rows);
            p_vals=zeros(c_rows);
            df=size(data,1);

            if c_cols ==1
                for i=1:c_rows 
                    o_col=data(:,col_1);
                    t_col=data(:,col_2);
                    [h,p]=ttest(o_col,t_col);

                    h_vals(i)=h;
                    p_vals(i)=p;
                end
                t_vals=tinv(p_vals,df);
            else
                disp('Colums out of range');
            end
        end
        
        %********************************************************************************
                 
        function WriteDataintoFile(data,outfile,outsepa)
            a=convertallfiletocell(data);
            dlmcell(outfile,a,outsepa);
        end
    end
    
    %######################********************************************************************************########################
    
    methods (Static,  Access = public)
        
        function ProbesetAveragedSignalTvalueOverExon(file,acols,sgnlcols,mergecols,ccols,Pval,fccutoff,log_base,shift,outfile)
            if(exist(file) == 0)
                disp(file,'does not exist');
            else
                a=importdata(file);
                data=convertallfiletocell(a);
                gene_sig=AveragedGeneSignalsOverTranscript(data,acols,sgnlcols,mergecols,log_base,shift);
                sort_data=sort_myown(data,acols(2),'ascending');
                disp(['data sorted']);
                signal = PartofCellString2Double(sort_data,sgnlcols);
                disp(['signal columns retrieved']);
                ids = PartofCellString2Double(sort_data,acols);
                row=size(data,1);
                row_m=size(mergecols,1);
                row_c= size(ccols,1);
                col_g=size(gene_sig,2);
                
                lcols(1:row_m)=1:row_m;
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
                                [~,p_val,t_val]=ma_analysis.TvalueCalculation(lv,ccols(L,1),ccols(L,2));
                            end

                            cnt=0;
                            for L=1:row_c
                                if(p_val <= Pval(L,1) && fc(L,1) > fccutoff)
                                    cnt=cnt+1;
                                end
                            end

                            if(cnt > 1)
                               sda(k,1)=ids(i,1);
                               sda(k,2)=ids(i,2);
                               for j=1:row_c
                                   sda(k,2+j)=fc(j,1);
                                   sda(k,2+row_c+j)=t_val(1,j);

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
                dlmwrite(outfile,sda,',');
            end
        end
        
        %********************************************************************************
        
        function [out] = AveragedGeneSignalsOverTranscript(filename,acols,sgnlcols,mergecols,log_base,shift)
            if(exist(filename) == 0)
                disp(filename,'does not exist');
            else
                data=importdata(filename);
                sort_data=sorting(data,acols,'ascending');
                signal=PartofCellString2Double(sort_data,sgnlcols);
                ids=PartofCellString2Double(sort_data,acols);
                [row,~] = size(data);
                row_cols=size(mergecols,1);
                exon_count=0; tv=[]; k=1; 
                lcols(1:row_cols)=1:row_cols;
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
            end
        end
    end
end