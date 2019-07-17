# Batch Driver for running buffalo simulations

NSIM=$1

for i in 2000 1800 1600 1400 1200 1000; do
  DATE=`date +%Y-%m-%d`
  nohup Rscript ../code/distance_sim_wrapper.R 0.5 ${i} ${NSIM} >> ${i}_dist_sim_${DATE}.log &
 done
