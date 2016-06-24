function [new_a] = sort_myown(indata,sort_colno,sort_type)


%indata=importdata('compare_trail2.txt'); sort_colno=5; sort_type='ascending';
[J]=convertallfiletocell(indata);

[row,col]=size(J);
s=sort_colno;


for i=1:row
    y=ischar(J{i,s});
    if(y==1)
        g=str2double(J{i,s});
        a(i,1)=g;
    else
        a(i,1)=(J{i,s}); 
    end
end
store=a;
[row_a,col_a]=size(a);

End=row_a;
Sim=ceil(row_a/2);
OddEven = mod(row_a,2); 

for j=1:Sim
    [MIN,min_i]=min(a(1,1)); [MAX,max_i]=max(a(1,1));
    [row_a,col_a]=size(a);

    start_a{j,:}=J{min_i,:};
    End_a{Sim-j+1,:}=J{max_i,:};

    K=1; A=[]; new_j=cell(1,1);
    for k=1:row_a
        if(k ~= min_i && k ~= max_i)
            A(K,1)=a(k,1);
            new_j{K,:}=J{k,:}; 
            K=K+1;
        end
    end
    a=zeros(size(A,1),size(A,2));
    J=zeros(size(A,1),size(A,2));
    a=A; J=new_j;
end

%ascending or descending

if(strcmp(sort_type,'ascending')==1)
    new_a=start_a; R=1;
    if(OddEven==0),Start=1; else Start=2; end
    for i=Start:size(End_a,1)
        new_a{Sim+R,:}=End_a{i,:};   
        R=R+1;
    end
    
else
    
    for i=1:size(End_a,1)
        for j=1:col
            new_a{i,j}=End_a{Sim-i+1,j};   
        end
    end
    if(OddEven==0),Start=Sim; else Start=Sim-1; end
    for i=1:Start
        new_a{Sim+i,:}=start_a{Start-i+1,:};   
    end
    
end



