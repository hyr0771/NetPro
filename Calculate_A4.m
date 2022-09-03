function [Drug_sim_corrf] = Calculate_A4(A,Drug_disease_data)
% A4

DR=Drug_disease_data';
RR=A;

[rows,cols] = size(DR);
y_d = zeros(rows, cols);            

col_no = find(sum(DR, 1) == 0);   
col_no=col_no';  

for i = 1 : length(col_no)
    RR(col_no(i), col_no(i)) = 0; 
    [sort_d, idx_d] = sort(RR(col_no(i), :), 'descend');  
    sum_d = sum(sort_d(1, 1 : 10));
    
    for j = 1 : 10          
        y_d(:, col_no(i)) =  y_d(:, col_no(i))+ sort_d(1, j) * DR(:,idx_d(1, j));
    end

    y_d(:, col_no(i)) = y_d(:,col_no(i)) / sum_d;
    
end
DR_new=DR+y_d;

%------------------------------------------------------------------------------------------------%
R=corrcoef(DR_new);
R_row = R(:)';  
R_row_normalizing = mapminmax(R_row, 0, 1); 
Mapped_R = reshape(R_row_normalizing, size(R)); 


%------------------------------------------------------------------------------------------------%

[rows, cols] = size( A );
knn_sim_drug = zeros(rows, cols);
[sort_network,idx]=sort(A,2,'descend');
K=10;
distan = zeros(rows, K);
[rows_t,cols_t]=find(DR_new'*DR_new);
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
     knn_sim_drug(i,j) =1/(max(distan(i,:))+max(distan(j,:)));
    end
end
knn_sim_drug(knn_sim_drug==inf)=0;


Drug_sim_corrf=3.0*knn_sim_drug.*abs(R);    
Drug_sim_corrf(isnan(Drug_sim_corrf))=0;
Drug_sim_corrf(logical(eye(size(Drug_sim_corrf))))=1;

