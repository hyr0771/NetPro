function [Disease_sim_corrf] = Calculate_S2(S1,Drug_disease_data)
% S2

RD=Drug_disease_data;
DD=S1;

[rows,cols] = size(RD);
y_d = zeros(rows, cols);              

col_no = find(sum(RD, 1) == 0);   
col_no=col_no'; 

for i = 1 : length(col_no)
    DD(col_no(i), col_no(i)) = 0;  
    [sort_d, idx_d] = sort(DD(col_no(i), :), 'descend');  
    sum_d = sum(sort_d(1, 1 : 10));
    
    for j = 1 : 10          
        y_d(:, col_no(i)) =  y_d(:, col_no(i))+ sort_d(1, j) * RD(:,idx_d(1, j));
    end
    
    y_d(:, col_no(i)) = y_d(:,col_no(i)) / sum_d;
    
end
RD_new=RD+y_d;



%------------------------------------------------------------------------------------------------%

R=corrcoef(RD_new);
R_row = R(:)';   
R_row_normalizing = mapminmax(R_row, 0, 1);  
Mapped_R = reshape(R_row_normalizing, size(R)); 


%------------------------------------------------------------------------------------------------%

[rows, cols] = size( S1 );
knn_sim_Disease = zeros(rows, cols);
[sort_network,idx]=sort(S1,2,'descend');
K=10;
distan = zeros(rows, K);
[rows_t,cols_t]=find(RD_new'*RD_new);
G = graph(rows_t,cols_t);
d = distances(G);
d(d==inf)=0;
for i=1:rows
    for j=1:K
        distan(i,j)=d(i,idx(i,j));
    end
end
for i=1:rows
    for j=1:rows
     knn_sim_Disease(i,j) =1/(max(distan(i,:))+max(distan(j,:)));
    end
end
knn_sim_Disease(knn_sim_Disease==inf)=0;


Disease_sim_corrf=3.0*knn_sim_Disease.*abs(R);    
Disease_sim_corrf(isnan(Disease_sim_corrf))=0;
Disease_sim_corrf(logical(eye(size(Disease_sim_corrf))))=1;

