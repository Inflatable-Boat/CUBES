#!/bin/sh
time -o times/v21_melt_12_p079_27 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  7.9 1.26610367 > logs/v21_melt_12_p079 &
time -o times/v21_compr_12_p079_27 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.40 750000 1000  7.9 1.26610367 > logs/v21_compr_12_p079 &
