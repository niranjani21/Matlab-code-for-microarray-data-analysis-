%function compare(filename_a,filename_b,colna_inarray,colnb_inarray)
%colna=colna_inarray;
clear;
clc;
colnb=[1];
[row_colnb,col_colnb]=size(colnb);
%colnb=colnb_inarray;
colna=[1];

[row_colna,col_colna]=size(colna);

a=importdata('compare_trail1.txt');
b=importdata('compare_trail2.txt');

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
for i=1:row_a
    for j=1:1
        y=ischar(a{i,j});
        if(y==1),g=str2double(a{i,j}); aval(i,j)=g; else aval(i,j)=(a{i,j}); end
    end
end
for i=1:row_b
    for j=colna(1):colnb(col_colnb)
        y=ischar(b{i,j});
        if(y==1),g=str2double(b{i,j}); bval(i,j)=g; else bval(i,j)=(b{i,j}); end
    end
end
% if row_a<row_b
%     row_c=row_b;
% else 
%     row_c=row_a;
% end
% col_c=col_a+col_b;
% c=cell(row_c,col_c);
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
    dlmcell('out_trail.txt',c,',');
end
