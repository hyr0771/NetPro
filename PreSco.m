function result_3=PreSco(AN, S, A5, S3, RD_adjmat_new)

lamta_1= 0.0001;
lamta_2= 0.05;

drpredict=lamta_1*(AN+A5)/5*RD_adjmat_new;
dipredict=lamta_2*RD_adjmat_new*(S+S3)/3;

result_3=0.5*drpredict+0.5*dipredict;
end