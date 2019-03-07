#!/bin/sh
/usr/bin/time -o times/v19_16_p010_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  1.0 1.26610367 > logs/v19_16_p010 &
/usr/bin/time -o times/v19_16_p030_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  3.0 1.26610367 > logs/v19_16_p030 &
/usr/bin/time -o times/v19_16_p050_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  5.0 1.26610367 > logs/v19_16_p050 &
/usr/bin/time -o times/v19_16_p070_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  7.0 1.26610367 > logs/v19_16_p070 &
/usr/bin/time -o times/v19_16_p080_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  8.0 1.26610367 > logs/v19_16_p080 &
/usr/bin/time -o times/v19_16_p100_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000 10.0 1.26610367 > logs/v19_16_p100 &
/usr/bin/time -o times/v19_16_p120_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000 12.0 1.26610367 > logs/v19_16_p120 &
/usr/bin/time -o times/v19_16_p150_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000 15.0 1.26610367 > logs/v19_16_p150 &

#/usr/bin/time -o times/v19_16_p071_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  7.1 1.26610367 > logs/v19_16_p071 &
/usr/bin/time -o times/v19_16_p072_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  7.2 1.26610367 > logs/v19_16_p072 &
#/usr/bin/time -o times/v19_16_p073_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  7.3 1.26610367 > logs/v19_16_p073 &
/usr/bin/time -o times/v19_16_p074_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  7.4 1.26610367 > logs/v19_16_p074 &
#/usr/bin/time -o times/v19_16_p075_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  7.5 1.26610367 > logs/v19_16_p075 &
/usr/bin/time -o times/v19_16_p076_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  7.6 1.26610367 > logs/v19_16_p076 &
#/usr/bin/time -o times/v19_16_p077_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  7.7 1.26610367 > logs/v19_16_p077 &
/usr/bin/time -o times/v19_16_p078_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  7.8 1.26610367 > logs/v19_16_p078 &
#/usr/bin/time -o times/v19_16_p079_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  7.9 1.26610367 > logs/v19_16_p079 &
/usr/bin/time -o times/v19_16_p090_27 -f "%E\t%U\t%S" ./v19_read_smartstart.exe c 16 0.50 500000 1000  9.0 1.26610367 > logs/v19_16_p090 &
