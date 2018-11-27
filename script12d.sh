#!/bin/sh
/usr/bin/time --output=lengths/time12d_p010_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 1.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time12d_p020_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 2.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time12d_p030_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 3.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time12d_p040_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 4.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time12d_p050_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 5.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time12d_p060_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 6.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time12d_p065_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 6.5  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time12d_p070_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 7.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time12d_p075_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 7.5  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time12d_p080_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 8.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time12d_p090_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 9.0  1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time12d_p100_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 10.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time12d_p120_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 12.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/time12d_p150_27.txt --format="%E real,\t%U user,\t%S sys" ./sl12d.exe 500001 0.50 15.0 1.26610367 > /dev/null &