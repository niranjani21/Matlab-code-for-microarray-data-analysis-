function [c]=sorting(data,sort_colno,sorttype)

[a]=convertallfiletocell(data);

if(strcmp('ascending',sorttype)==1)
    s=sort_colno;
else
    S=sprintf('-%d',sort_colno);
    s=str2double(S);
end

c=sortrows(a,s);
