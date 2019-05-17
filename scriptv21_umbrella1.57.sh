#!/bin/sh
time -o times/v30_20p6.17a1.57 -f "%E\t%U\t%S,\t%P\n" ./v30_sim.exe r datafolder/v32_20pf0.30p01.0a1.57 v21_20p6.17a1.57 750000 1000 6.17 1 > logs/v21_20p6.17a1.57 &
time -o times/v30_20p6.18a1.57 -f "%E\t%U\t%S,\t%P\n" ./v30_sim.exe r datafolder/v32_20pf0.30p01.0a1.57 v21_20p6.18a1.57 750000 1000 6.18 1 > logs/v21_20p6.18a1.57 &
time -o times/v30_20p6.19a1.57 -f "%E\t%U\t%S,\t%P\n" ./v30_sim.exe r datafolder/v32_20pf0.30p01.0a1.57 v21_20p6.19a1.57 750000 1000 6.19 1 > logs/v21_20p6.19a1.57 &
time -o times/v30_20p6.20a1.57 -f "%E\t%U\t%S,\t%P\n" ./v30_sim.exe r datafolder/v32_20pf0.30p01.0a1.57 v21_20p6.20a1.57 750000 1000 6.20 1 > logs/v21_20p6.20a1.57 &
time -o times/v30_20p6.25a1.57 -f "%E\t%U\t%S,\t%P\n" ./v30_sim.exe r datafolder/v32_20pf0.30p01.0a1.57 v21_20p6.20a1.57 750000 1000 6.20 1 > logs/v21_20p6.20a1.57 &
