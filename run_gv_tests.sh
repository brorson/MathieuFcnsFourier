#!/bin/bash

#  fns=("mathieu_ce_gvs_q-1.csv" "mathieu_ce_gvs_q-0.1.csv" "mathieu_ce_gvs_q0.1.csv" "mathieu_ce_gvs_q1.csv")

#------------------------------------------------
# Mathieu ce
fns=( $(ls *.csv | grep ce_gvs) )

for fn in ${fns[@]}
do
echo "Running tests using filename" "$fn"
echo "$fn" | octave test_mathieu_ce_gvs.m
done
echo ==========================================


#------------------------------------------------
# Mathieu se
fns=( $(ls *.csv | grep se_gvs) )

for fn in ${fns[@]}
do
echo "Running tests using filename" "$fn"
echo "$fn" | octave test_mathieu_se_gvs.m
done
echo ==========================================


#------------------------------------------------
# Mathieu ce deriv
fns=( $(ls *.csv | grep ce_deriv_gvs) )

for fn in ${fns[@]}
do
echo "Running tests using filename" "$fn"
echo "$fn" | octave test_mathieu_ce_deriv_gvs.m
done
echo ==========================================

#------------------------------------------------
# Mathieu se deriv
fns=( $(ls *.csv | grep se_deriv_gvs) )

for fn in ${fns[@]}
do
echo "Running tests using filename" "$fn"
echo "$fn" | octave test_mathieu_se_deriv_gvs.m
done
echo ==========================================

