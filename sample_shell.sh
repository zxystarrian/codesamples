################################################################################
# Generate PBS job scripts with varied combinations of parameters:
# 	path: No. of pathway
# 	part: Noise in the simulation, increase by 0.25
# 	thres: Threshold used in the EMVC R package
# 	cohorts
################################################################################

PBS_O_WORKDIR=/dartfs-hpc/rc/lab/F/FrostH/members/xzheng/project4/tumor/simulation/sim_partcor/nocluster_univar/cor0.2/path27

for path in {27..27}
do
for part in {0..0}
do
for thres in {"0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","0.91","0.93","0.95","0.97","0.99"}
do
for cohorts in {"BLCA","BRCA","CESC","COAD","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","PRAD","READ","SKCM","STAD","TGCT","THCA","UCEC"}
do
now=$(echo "$part*0.25"|bc)
touch $cohorts"_P_expr_"$thres$now$path.pbs
echo "#!/bin/bash -l" >> $cohorts"_P_expr_"$thres$now$path.pbs
echo "#PBS -l walltime=24:00:00" >> $cohorts"_P_expr_"$thres$now$path.pbs
echo "cd $PBS_O_WORKDIR" >> $cohorts"_P_expr_"$thres$now$path.pbs
echo "Rscript sim_optimized.R $cohorts $thres $now $path > P_expr_$thres$cohorts$now$path.out" >> $cohorts"_P_expr_"$thres$now$path.pbs
done
done
done
done
