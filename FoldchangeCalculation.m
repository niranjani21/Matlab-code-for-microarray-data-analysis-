function [fc] = FoldchangeCalculation(data,col_1,col_2)
%a=importdata('out.txt'); col_1=7; col_2=8;
%data=convertallfiletocell(a);

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