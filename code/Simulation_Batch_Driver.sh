# Batch Driver for running buffalo simulations

NSIM=$1

for i in {1..10}; do
  RHO=$(echo "scale = 1; ${i}/10" | bc);
  DATE=`date +%Y-%m-%d`
  nohup Rscript ../code/simulate_wrapper.R ${RHO} 10 ${NSIM} >> ${i}_sim_${DATE}.log &
 done
