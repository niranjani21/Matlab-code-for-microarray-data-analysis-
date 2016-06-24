function [t_val] = TvalueCalculation(data,col_1,col_2)

o_col=data(:,col_1);
t_col=data(:,col_2);

o_mean=mean(o_col);
t_mean=mean(t_col);
o_var=var(o_col);
t_var=var(t_col);
if(o_var >0 || t_var > 0)
    t_val = (o_mean-t_mean)/sqrt((o_var+t_var)/size(data,1));
end