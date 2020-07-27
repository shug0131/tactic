#!/bin/bash

scenario=$1
max=32
parallel=$2 
parallel=$(( $parallel < $max ? $parallel : $max ))
echo "no of cores is $parallel"
mkdir "/home/sjb277/rds/hpc-work/$scenario"
mkdir "/home/sjb277/rds/hpc-work/$scenario/output_files"
current_time=$(date "+%Y.%m.%d-%H.%M.%S")
cp /home/sjb277/tactic/future_real_data.R /home/sjb277/rds/hpc-work/$scenario/r_code_future_real_data$current_time.R
sbatch --export=scenario=$scenario,parallel=$parallel --output /home/sjb277/rds/hpc-work/$scenario/output_files/slurm-%A_%a.out --cpus-per-task=$parallel slurm_future_inference_skylake.txt