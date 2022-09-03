function RD_mat_new = RDKNN( RD_mat, RR_mat, DD_mat, K, r, t1, t2)     % RDKNN

[rows,cols]=size(RD_mat);          
y_m=zeros(rows,cols);  
y_d=zeros(rows,cols);  

%========================================================================================%

knn_network_m = KNN( RR_mat, K );       
for i = 1 : rows   
         w=zeros(1,K);                  
        [sort_m,idx_m]=sort(knn_network_m(i,:),2,'descend'); 
        sum_m=sum(sort_m(1,1:K));   
        for j = 1 : K
            w(1,j)=r^(j-1)*sort_m(1,j); 
            y_m(i,:) =  y_m(i,:)+ w(1,j)* RD_mat(idx_m(1,j),:); 
        end                      
            y_m(i,:)=y_m(i,:)/sum_m;              
end

%===================================================================================%

knn_network_d = KNN( DD_mat , K );       
for i = 1 : cols   
        w=zeros(1,K);
        [sort_d,idx_d]=sort(knn_network_d(i,:),2,'descend');
        sum_d=sum(sort_d(1,1:K));
        for j = 1 : K
            w(1,j)=r^(j-1)*sort_d(1,j);
            y_d(:,i) =  y_d(:,i)+ w(1,j)* RD_mat(:,idx_d(1,j)); 
        end                      
            y_d(:,i)=y_d(:,i)/sum_d;               
end

%====================================================================================%

y_md=y_m*t1+y_d*t2;  
RD_mat_new=max(RD_mat,y_md);

end

%========================================================================================%

function [ knn_network ] = KNN( network , k )
    [rows, cols] = size( network );
    network= network-diag(diag(network));    
    knn_network = zeros(rows, cols);
    [sort_network,idx]=sort(network,2,'descend'); 
    
    for i = 1 : rows
        knn_network(i,idx(i,1:k))=sort_network(i,1:k);
     
    end
    
end


