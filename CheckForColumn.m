function [cnt] = CheckForColumn(col_no,check_cols)
%col_no=7; check_cols=sig;
col=size(check_cols,2);
for i=1:col
    if(check_cols(i)==col_no), cnt(1)=1; cnt(2)=i; break; else cnt(1)=0; cnt(2)=0; end
end
