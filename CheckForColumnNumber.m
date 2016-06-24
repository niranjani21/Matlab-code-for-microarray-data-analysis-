function [cnt,sc] = CheckForColumnNumber(col_no,check_cols)
%col_no=7; check_cols=sig;
col=size(check_cols,2);
for i=1:col
    if(check_cols(i)==col_no),sc=i; cnt=1; break; else cnt=0;end
end