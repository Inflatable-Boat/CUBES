#!/bin/sh
/usr/bin/time -o lengths/10p010_27.txt -f "%E real,\t%U user,\t%S sys" ./v1read.exe c 10 500001 0.50 1.0  1.26610367 > logs/v1_10_p010.log &
/usr/bin/time -o lengths/10p030_27.txt -f "%E real,\t%U user,\t%S sys" ./v1read.exe c 10 500001 0.50 3.0  1.26610367 > logs/v1_10_p030.log &
/usr/bin/time -o lengths/10p050_27.txt -f "%E real,\t%U user,\t%S sys" ./v1read.exe c 10 500001 0.50 5.0  1.26610367 > logs/v1_10_p050.log &
/usr/bin/time -o lengths/10p070_27.txt -f "%E real,\t%U user,\t%S sys" ./v1read.exe c 10 500001 0.50 7.0  1.26610367 > logs/v1_10_p070.log &
/usr/bin/time -o lengths/10p080_27.txt -f "%E real,\t%U user,\t%S sys" ./v1read.exe c 10 500001 0.50 8.0  1.26610367 > logs/v1_10_p080.log &
/usr/bin/time -o lengths/10p100_27.txt -f "%E real,\t%U user,\t%S sys" ./v1read.exe c 10 500001 0.50 10.0 1.26610367 > logs/v1_10_p100.log &
/usr/bin/time -o lengths/10p120_27.txt -f "%E real,\t%U user,\t%S sys" ./v1read.exe c 10 500001 0.50 12.0 1.26610367 > logs/v1_10_p120.log &
/usr/bin/time -o lengths/10p150_27.txt -f "%E real,\t%U user,\t%S sys" ./v1read.exe c 10 500001 0.50 15.0 1.26610367 > logs/v1_10_p150.log &
