#!/bin/sh

/usr/bin/time -o times/v21_melt_12_p030_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  3.0 1.57079632679 > logs/v21_melt_12_p030_57 &
/usr/bin/time -o times/v21_melt_12_p040_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  4.0 1.57079632679 > logs/v21_melt_12_p040_57 &
/usr/bin/time -o times/v21_melt_12_p045_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  4.5 1.57079632679 > logs/v21_melt_12_p045_57 &
/usr/bin/time -o times/v21_melt_12_p050_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  5.0 1.57079632679 > logs/v21_melt_12_p050_57 &

/usr/bin/time -o times/v21_melt_12_p052_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  5.2 1.57079632679 > logs/v21_melt_12_p052_57 &
/usr/bin/time -o times/v21_melt_12_p054_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  5.4 1.57079632679 > logs/v21_melt_12_p054_57 &
/usr/bin/time -o times/v21_melt_12_p056_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  5.6 1.57079632679 > logs/v21_melt_12_p056_57 &
/usr/bin/time -o times/v21_melt_12_p058_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  5.8 1.57079632679 > logs/v21_melt_12_p058_57 &

/usr/bin/time -o times/v21_melt_12_p060_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  6.0 1.57079632679 > logs/v21_melt_12_p060_57 &
/usr/bin/time -o times/v21_melt_12_p065_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  6.5 1.57079632679 > logs/v21_melt_12_p065_57 &

/usr/bin/time -o times/v21_melt_12_p070_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  7.0 1.57079632679 > logs/v21_melt_12_p070_57 &
/usr/bin/time -o times/v21_melt_12_p080_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000  8.0 1.57079632679 > logs/v21_melt_12_p080_57 &
/usr/bin/time -o times/v21_melt_12_p100_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000 10.0 1.57079632679 > logs/v21_melt_12_p100_57 &
/usr/bin/time -o times/v21_melt_12_p120_57 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe c 12 0.60 750000 1000 12.0 1.57079632679 > logs/v21_melt_12_p120_57 &
