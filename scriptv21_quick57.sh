#!/bin/sh
/usr/bin/time -o times/v21_compr_12_p061_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.40 750000 1000  6.1 1.57079632679 > logs/v21_compr_12_p061_57 &
/usr/bin/time -o times/v21_compr_12_p062_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.40 750000 1000  6.2 1.57079632679 > logs/v21_compr_12_p062_57 &
/usr/bin/time -o times/v21_compr_12_p063_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.40 750000 1000  6.3 1.57079632679 > logs/v21_compr_12_p063_57 &
/usr/bin/time -o times/v21_compr_12_p064_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.40 750000 1000  6.4 1.57079632679 > logs/v21_compr_12_p064_57 &
/usr/bin/time -o times/v21_melt_12_p061_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  6.1 1.57079632679 > logs/v21_melt_12_p061_57 &
/usr/bin/time -o times/v21_melt_12_p062_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  6.2 1.57079632679 > logs/v21_melt_12_p062_57 &
/usr/bin/time -o times/v21_melt_12_p063_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  6.3 1.57079632679 > logs/v21_melt_12_p063_57 &
/usr/bin/time -o times/v21_melt_12_p064_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  6.4 1.57079632679 > logs/v21_melt_12_p064_57 &
