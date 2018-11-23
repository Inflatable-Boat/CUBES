#!/bin/sh
/usr/bin/time --output=lengths/time10p010_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10_cl.exe 200001 0.50 1.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10p030_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10_cl.exe 200001 0.50 3.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10p050_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10_cl.exe 200001 0.50 5.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10p070_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10_cl.exe 200001 0.50 7.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10p080_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10_cl.exe 200001 0.50 8.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10p100_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10_cl.exe 200001 0.50 10.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10p120_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10_cl.exe 200001 0.50 12.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10p150_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10_cl.exe 200001 0.50 15.0 1.26610367 > /dev/null &
