#!/bin/sh
/usr/bin/time --output=lengths/time10d_6fix_p010_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 1.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10d_6fix_p020_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 2.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10d_6fix_p030_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 3.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10d_6fix_p040_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 4.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10d_6fix_p050_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 5.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10d_6fix_p060_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 6.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10d_6fix_p065_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 6.5  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10d_6fix_p070_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 7.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10d_6fix_p075_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 7.5  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10d_6fix_p080_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 8.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10d_6fix_p090_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 9.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10d_6fix_p100_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 10.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10d_6fix_p120_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 12.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time10d_6fix_p150_27.txt --format="%E real,\t%U user,\t%S sys" ./sl10d_6fix.exe 500001 0.50 15.0 1.26610367 > /dev/null &