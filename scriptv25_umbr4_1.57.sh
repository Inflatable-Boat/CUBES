#!/bin/sh
time -o times/v25_20p6.17a1.57_step3_t200 -f "%E\t%U\t%S,\t%P\n" ./v25_sim_umbrella.exe r datafolder/v25_20p6.17a1.57_step2_t100/coords_step0207000.poly v25_20p6.17a1.57_step3_t200 300000 1000 6.17 1 t 200 0.005 > logs/v25_20p6.17a1.57_step3_t200 &
time -o times/v25_20p6.17a1.57_step3_t250 -f "%E\t%U\t%S,\t%P\n" ./v25_sim_umbrella.exe r datafolder/v25_20p6.17a1.57_step2_t150/coords_step0284000.poly v25_20p6.17a1.57_step3_t250 300000 1000 6.17 1 t 250 0.005 > logs/v25_20p6.17a1.57_step3_t250 &
time -o times/v25_20p6.17a1.57_step3_t300 -f "%E\t%U\t%S,\t%P\n" ./v25_sim_umbrella.exe r datafolder/v25_20p6.17a1.57_step2_t200/coords_step0266000.poly v25_20p6.17a1.57_step3_t300 300000 1000 6.17 1 t 300 0.005 > logs/v25_20p6.17a1.57_step3_t300 &
