#!/bin/sh
/usr/bin/time -o lengths/time10p010_27.txt -f "%E real,\t%U user,\t%S sys" ./sl10_cl.exe 500001 0.50 1.0  1.26610367 > logs/sl10_p010.log &
/usr/bin/time -o lengths/time10p030_27.txt -f "%E real,\t%U user,\t%S sys" ./sl10_cl.exe 500001 0.50 3.0  1.26610367 > logs/sl10_p030.log &
/usr/bin/time -o lengths/time10p050_27.txt -f "%E real,\t%U user,\t%S sys" ./sl10_cl.exe 500001 0.50 5.0  1.26610367 > logs/sl10_p050.log &
/usr/bin/time -o lengths/time10p070_27.txt -f "%E real,\t%U user,\t%S sys" ./sl10_cl.exe 500001 0.50 7.0  1.26610367 > logs/sl10_p070.log &
/usr/bin/time -o lengths/time10p080_27.txt -f "%E real,\t%U user,\t%S sys" ./sl10_cl.exe 500001 0.50 8.0  1.26610367 > logs/sl10_p080.log &
/usr/bin/time -o lengths/time10p100_27.txt -f "%E real,\t%U user,\t%S sys" ./sl10_cl.exe 500001 0.50 10.0 1.26610367 > logs/sl10_p100.log &
/usr/bin/time -o lengths/time10p120_27.txt -f "%E real,\t%U user,\t%S sys" ./sl10_cl.exe 500001 0.50 12.0 1.26610367 > logs/sl10_p120.log &
/usr/bin/time -o lengths/time10p150_27.txt -f "%E real,\t%U user,\t%S sys" ./sl10_cl.exe 500001 0.50 15.0 1.26610367 > logs/sl10_p150.log &
