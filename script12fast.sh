#!/bin/sh
/usr/bin/time -o times/12p010_27 -f "%E\t%U\t%S" ./v1read c 12 500001 0.50  1.0 1.26610367 > logs/v1_12_p010 &
/usr/bin/time -o times/12p030_27 -f "%E\t%U\t%S" ./v1read c 12 500001 0.50  3.0 1.26610367 > logs/v1_12_p030 &
/usr/bin/time -o times/12p050_27 -f "%E\t%U\t%S" ./v1read c 12 500001 0.50  5.0 1.26610367 > logs/v1_12_p050 &
/usr/bin/time -o times/12p070_27 -f "%E\t%U\t%S" ./v1read c 12 500001 0.50  7.0 1.26610367 > logs/v1_12_p070 &
/usr/bin/time -o times/12p080_27 -f "%E\t%U\t%S" ./v1read c 12 500001 0.50  8.0 1.26610367 > logs/v1_12_p080 &
/usr/bin/time -o times/12p100_27 -f "%E\t%U\t%S" ./v1read c 12 500001 0.50 10.0 1.26610367 > logs/v1_12_p100 &
/usr/bin/time -o times/12p120_27 -f "%E\t%U\t%S" ./v1read c 12 500001 0.50 12.0 1.26610367 > logs/v1_12_p120 &
/usr/bin/time -o times/12p150_27 -f "%E\t%U\t%S" ./v1read c 12 500001 0.50 15.0 1.26610367 > logs/v1_12_p150 &
