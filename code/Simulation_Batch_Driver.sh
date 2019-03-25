# Batch Driver for running buffalo simulations

for i in {1..10}; do
  RHO=$(echo "scale = 1; ${i}/10" | bc);
  DATE=`date +%Y-%m-%d`
  nohup RScript ../code/simulate_wrapper.R ${RHO} 10 10 >> ${i}_sim_${DATE}.log &
 done
