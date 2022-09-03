clear
clc
warning('off');
addpath('Datasets');

%================================================================================%
% The PREDICT dataset

[Disease_drug]=importdata('predict_admat_dgc.txt'); 
Disease_drug_data=Disease_drug.data;                   
Drug_disease_data=Disease_drug_data';

[Disease_sim]=importdata('predict_simmat_dg.txt'); 
Disease_sim_data=Disease_sim.data;                    
S1=Disease_sim_data;                                    
                         
[Input]=importdata('predict_simmat_dc_go.txt'); 
drug_go=Input.data;
A1=drug_go;
[Input]=importdata('predict_simmat_dc_domain.txt'); 
drug_domain=Input.data;
A2=drug_domain;
[Input]=importdata('predict_simmat_dc_chemical.txt'); 
drug_chemical=Input.data;
A3=drug_chemical;

%====================================================================================%
% The fourth drug similarity and the second disease similarity are calculated here

[S2] = Calculate_S2(S1,Drug_disease_data); 
[Drug_sim_A1_corrf] = Calculate_A4(A1,Drug_disease_data);
[Drug_sim_A2_corrf] = Calculate_A4(A2,Drug_disease_data);
[Drug_sim_A3_corrf] = Calculate_A4(A3,Drug_disease_data);
A4=(Drug_sim_A1_corrf+Drug_sim_A2_corrf+Drug_sim_A3_corrf)/3;

AN=A1+A2+A3+A4;
S=S1+S2;

%==================================================================================%

KK=10;
r=0.5;
t1=0.5;
t2=0.5;
RD_adjmat_new=RDKNN( Drug_disease_data, (AN/4), (S/2), KK, r, t1, t2);

interaction_matrix=RD_adjmat_new'; 
neighbor_num=100;
alpha=0.5;
similairty_matrix_disease=Calculate_Similairty(interaction_matrix,0,neighbor_num,'regulation2'); 
LP_disease_result=Label_Propagation(similairty_matrix_disease,Drug_disease_data',alpha);         
LP_disease_result=LP_disease_result';  % from the disease side


interaction_matrix=RD_adjmat_new;   
neighbor_num=200;
alpha=0.5;
similairty_matrix_drug=Calculate_Similairty(interaction_matrix,0,neighbor_num,'regulation2');    
LP_Drug_result=Label_Propagation(similairty_matrix_drug,Drug_disease_data,alpha); % from the drug side

w=0.7;
result_1=(1-w)*(LP_disease_result)+w*LP_Drug_result;

%=====================================================================================%

KK=10;
r=1.0;
t1=0.10;
t2=0.15;
result_2=RDKNN( Drug_disease_data, (AN/4), (S/2), KK, r, t1, t2);

%====================================================================================%

[A5,S3]=NewSim_A5_S3(Drug_disease_data);
result_3=PreSco(AN, S, A5, S3, RD_adjmat_new);

%==================================================================================%
% Integrating Information And Normalizing

M_recovery=result_1+result_2+result_3;
M_recovery(M_recovery<0)=0;
M_recovery(M_recovery>1)=1;



