#!/bin/sh
/usr/bin/time -o times/10p010_57 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50  1.0 1.57079632679 > logs/v1_10_p010_57 &
/usr/bin/time -o times/10p030_57 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50  3.0 1.57079632679 > logs/v1_10_p030_57 &
/usr/bin/time -o times/10p050_57 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50  5.0 1.57079632679 > logs/v1_10_p050_57 &
/usr/bin/time -o times/10p070_57 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50  7.0 1.57079632679 > logs/v1_10_p070_57 &
/usr/bin/time -o times/10p080_57 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50  8.0 1.57079632679 > logs/v1_10_p080_57 &
/usr/bin/time -o times/10p100_57 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50 10.0 1.57079632679 > logs/v1_10_p100_57 &
/usr/bin/time -o times/10p120_57 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50 12.0 1.57079632679 > logs/v1_10_p120_57 &
/usr/bin/time -o times/10p150_57 -f "%E\t%U\t%S" ./v1read c 10 500001 0.50 15.0 1.57079632679 > logs/v1_10_p150_57 &
