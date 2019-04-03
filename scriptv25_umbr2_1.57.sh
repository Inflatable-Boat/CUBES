#!/bin/sh
time -o times/v25_20p6.17a1.57_step2_t050 -f "%E\t%U\t%S,\t%P\n" ./v25_sim_umbrella.exe r start_positions/v21_20p6.17a1.57.poly v25_20p6.17a1.57_step2_t050 300000 1000 6.17 1 t  50 0.005 > logs/v25_20p6.17a1.57_step2_t050 &
time -o times/v25_20p6.17a1.57_step2_t200 -f "%E\t%U\t%S,\t%P\n" ./v25_sim_umbrella.exe r start_positions/v21_20p6.17a1.57.poly v25_20p6.17a1.57_step2_t200 300000 1000 6.17 1 t 200 0.005 > logs/v25_20p6.17a1.57_step2_t200 &
