# graph type: band,er,knn,sf
g = knn 
# model; it depends on how $\lambda$ is generated
m = glasso 
# random seed
nt = 1 
# method for tuning: K#CV, LOOCV, ACV, ACV_Kfold, ProxACV, BIC, EBIC
method = K3CV 
# the number of samples
n = 100 
# the number of nodes
p = 10 

# Generating data; need to set global variables N_VAL, P_VAL, G_VAL, DIST_TYPE
Rscript gen_data.R ${nt}

# Generating $\lambda$; default number of $\lambda$ is 100
Rscript gen_lam.R --g={g} --nt={nt}

# Testing a specific method; results saved to acv_save/
Rscript main.R --n={n} --p={p} --g={g} --m={m} --nt={nt} --method={method}

# Process post-cleaning; results saved to EXP_FDR
Rscript postcleaning.R --n={n} --p={p} --g={g} --m={m} --nt={nt}

