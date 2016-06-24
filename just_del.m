a=importdata('abc.txt');
b=sort(a(:,1));
disp(b);



% tic
% z=importdata('HG-U133_probesetid_x_y.txt');
% 
% for d=1:length(z)-1
%         a=z(d,2);
%         
%         for j=d+1:length(z)
%             if a<z(j,2)
%                 z(d,2) = a;
%                 
%                 
%             else
%                 a=z(j,2);
%                 a1=z(j,3);
%                 a2=z(j,1);
%                 b=z(d,2); 
%                 b1=z(d,3);
%                 b2=z(d,1);
%                 z(d,2)=a;
%                 z(d,3)=a1;
%                 z(d,1)=a2;
%                 z(j,2)=b;
%                 z(j,3)=b1;
%                 z(j,1)=b2;
%             end
%         end
% end
% 
% fsav=fopen('sort_probeset.txt','wt');
% for j=1:length(z)
%     fprintf(fsav,'%d,%d,%d\n',z(j,1),z(j,2),z(j,3));
% end
% fclose(fsav);

% 
% a=importdata('comb_cere_heart_kid.txt');
% b=importdata('sort_probeset.txt');
% c=[];
% j=1;
% for i=1:length(b)
%     while j<length(b)
%         if (a(i,1) < b(j,2))
%             break;
%         else if (a(i,1) == b(j,2))
%                 c(i,1)=a(i,1);
%                 c(i,2)=a(i,2);
%                 c(i,3)=a(i,3);
%                 c(i,4)=b(j,1);
%                 j=j+1;
%             else if (a(i,1) > b(j,1))
%                     j=j+1;
%                 end
%             end
%         end
%     end
% end
% 
% fsav=fopen('commondata.txt','wt');
% for i=1:length(c)
%     fprintf(fsav,'%d,%d,%d,%d\n',c(i,1),c(i,2),c(i,3),c(i,4));
% end
% fclose(fsav);

% for i=1:100
%     for j=1:100
%         if (a.data(i,1)< b.data(i,2))
%             break;
%         else if (a(i,1) == b(i,2))
%                 c(i,1)=a(i,3);
%                 c(i,2)=b(i,1);
%             
%             end
%         end
%     end
% end
% 
% fsav = fopen('out.txt', 'wt');
% fprintf(fsav,'c(i,1)c(i,2)\n');
% for i=1:2
%     fprintf(fsav,'%d %d\n',c(i,1),c(i,2));
% end
% fclose(fsav);
% 
% 
% % (Vectormnk<char> *v, int left, int right, int col_no, int comp_stat) {
% 	i=left;  j=right;
% 	Vectorm<char> tempx_vect(v->k), *tx = &tempx_vect;
%     
% 	for k=0:length(v)
% 		tx(k,1) = v((i+j)/2,  v->x[static_cast<int> ((i + j)/2)][col_no][k];
% 	do { 
% 		if(comp_stat == 0) {
% 			while(atof(v->x[i][col_no]) < atof(tx->x))
% 				i++; 
% 			while(atof(v->x[j][col_no]) > atof(tx->x))
% 				j--; 
% 		}
% 		else {
% 			while(strcmp(v->x[i][col_no], tx->x) < 0)
% 				i++; 
% 			while(strcmp(v->x[j][col_no], tx->x) > 0)
% 				j--;
% 		}
% 		if (i <= j) { 
% 			swapvectormnkrows(v,i,j); 
% 			i++; 
% 			j--;  
% 		}
% 	}
% 	while (i <= j);
% 	if(left < j) 
% 		QuickSortVectormnkWRTNColumnAscend(v, left, j, col_no,comp_stat);
% 	if(i < right) 
% 		QuickSortVectormnkWRTNColumnAscend(v, i, right, col_no,comp_stat);