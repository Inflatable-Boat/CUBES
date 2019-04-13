#!/bin/sh
time -o times/v25_20p6.17a1.57_step1_t050long -f "%E\t%U\t%S,\t%P\n" ./c25_sim_umbrella_noflush.exe r umbr_starts/start1.poly v25_20p6.17a1.57_step1_t050long 1000000 1000 6.17 1 t 50 0.005 > logs/v25_20p6.17a1.57_step1_t050long &
time -o times/v25_20p6.17a1.57_step2_t070long -f "%E\t%U\t%S,\t%P\n" ./c25_sim_umbrella_noflush.exe r umbr_starts/start1.poly v25_20p6.17a1.57_step2_t070long 1000000 1000 6.17 1 t 70 0.005 > logs/v25_20p6.17a1.57_step2_t070long &
time -o times/v25_20p6.17a1.57_step3_t090long -f "%E\t%U\t%S,\t%P\n" ./c25_sim_umbrella_noflush.exe r umbr_starts/start1.poly v25_20p6.17a1.57_step3_t090long 1000000 1000 6.17 1 t 90 0.005 > logs/v25_20p6.17a1.57_step3_t090long &
