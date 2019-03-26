#!/bin/sh
time -o times/v21_20p7.60a1.27 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.27.poly v21_20p7.60a1.27 750000 1000 7.60 1 > logs/v21_20p7.60a1.27 &
time -o times/v21_20p7.65a1.27 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.27.poly v21_20p7.65a1.27 750000 1000 7.65 1 > logs/v21_20p7.65a1.27 &
time -o times/v21_20p7.70a1.27 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.27.poly v21_20p7.70a1.27 750000 1000 7.70 1 > logs/v21_20p7.70a1.27 &
time -o times/v21_20p7.75a1.27 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.27.poly v21_20p7.75a1.27 750000 1000 7.75 1 > logs/v21_20p7.75a1.27 &
time -o times/v21_20p7.80a1.27 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.27.poly v21_20p7.80a1.27 750000 1000 7.80 1 > logs/v21_20p7.80a1.27 &
time -o times/v21_20p7.85a1.27 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.27.poly v21_20p7.85a1.27 750000 1000 7.85 1 > logs/v21_20p7.85a1.27 &
time -o times/v21_20p7.90a1.27 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.27.poly v21_20p7.90a1.27 750000 1000 7.90 1 > logs/v21_20p7.90a1.27 &
time -o times/v21_20p7.95a1.27 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.27.poly v21_20p7.95a1.27 750000 1000 7.95 1 > logs/v21_20p7.95a1.27 &
