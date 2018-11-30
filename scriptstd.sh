#!/bin/sh
/usr/bin/time -o times/10p010_27 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50  1.0 1.26610367 > logs/v1_10_p010 &
/usr/bin/time -o times/10p030_27 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50  3.0 1.26610367 > logs/v1_10_p030 &
/usr/bin/time -o times/10p050_27 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50  5.0 1.26610367 > logs/v1_10_p050 &
/usr/bin/time -o times/10p070_27 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50  7.0 1.26610367 > logs/v1_10_p070 &
/usr/bin/time -o times/10p080_27 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50  8.0 1.26610367 > logs/v1_10_p080 &
/usr/bin/time -o times/10p100_27 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50 10.0 1.26610367 > logs/v1_10_p100 &
/usr/bin/time -o times/10p120_27 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50 12.0 1.26610367 > logs/v1_10_p120 &
/usr/bin/time -o times/10p150_27 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50 15.0 1.26610367 > logs/v1_10_p150 &
