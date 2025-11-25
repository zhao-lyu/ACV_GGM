# graph type: band,er,knn,sf
g=band
# model; it depends on how $\lambda$ is generated
m=glasso
# random seed
nt=1
# method for tuning: K#CV, LOOCV, ACV, ACV_Kfold, ProxACV, BIC, EBIC
# required for demo: K3CV, BIC, ACV_Kfold
methods=(K3CV BIC ACV_Kfold)
# the number of samples
n=100
# the number of nodes
p=10

dist_type=lognormal

# Generating data; need to set global variables N_VAL, P_VAL, G_VAL, DIST_TYPE
# set global variables N_VAL, P_VAL, G_VAL, DIST_TYPE
export N_VAL=${n}
export P_VAL=${p}
export G_VAL=${g}
export DIST_TYPE=${dist_type}
Rscript gen_data.R ${nt}

# Generating $\lambda$; default number of $\lambda$ is 100
Rscript gen_lam.R --g=${g} --nt=${nt}

# Testing methods in loop; results saved to acv_save/
for method in "${methods[@]}"; do
    echo "Running method: ${method}"
    Rscript main.R --n=${n} --p=${p} --g=${g} --m=${m} --nt=${nt} --method=${method}
done

# Process post-cleaning; results saved to EXP_FDR
Rscript postcleaning.R --n=${n} --p=${p} --g=${g} --m=${m} --nt=${nt}

