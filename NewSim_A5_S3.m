function [Wrr_new,Wdd_new]=NewSim_A5_S3(Drug_disease_data)

    
[row_num,col_num]=size(Drug_disease_data);

Wrr=zeros(row_num,row_num);
for i=1:row_num
        
        for j=i:row_num
            Wrr(i,j)=sum((Drug_disease_data(i,:)+Drug_disease_data(j,:)==2));
            Wrr(j,i)=Wrr(i,j);
        end
               
        
end
    
Wrr_new=zeros(row_num,row_num);
for i=1:row_num      
           
    Wrr_new(i,:)=Wrr(i,:)/col_num;
                   
end
Wrr_new(find(isnan(Wrr_new)==1))=0;

Wrr_new=Wrr_new*70;       % A5


%-----------------------------------------------------------------------------%

Wdd=zeros(col_num,col_num);
for i=1:col_num
        
        for j=i:col_num
            Wdd(i,j)=sum((Drug_disease_data(:,i)+Drug_disease_data(:,j)==2));
            Wdd(j,i)=Wdd(i,j);
        end
               
        
end

Wdd_new=zeros(col_num,col_num);
for i=1:col_num       
           
    Wdd_new(i,:)=Wdd(i,:)/row_num;
                 
end
Wdd_new(find(isnan(Wdd_new)==1))=0;

Wdd_new=Wdd_new*70;     % S3

end


