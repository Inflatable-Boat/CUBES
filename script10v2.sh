#!/bin/sh
/usr/bin/time --output=lengths/time10clv2_p010_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10cl.exe 500001 0.50 1.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10clv2_p030_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10cl.exe 500001 0.50 3.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10clv2_p050_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10cl.exe 500001 0.50 5.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10clv2_p070_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10cl.exe 500001 0.50 7.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10clv2_p080_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10cl.exe 500001 0.50 8.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10clv2_p100_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10cl.exe 500001 0.50 10.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10clv2_p120_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10cl.exe 500001 0.50 12.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10clv2_p150_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10cl.exe 500001 0.50 15.0 1.26610367 > /dev/null &
