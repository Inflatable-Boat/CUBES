#!/bin/sh
/usr/bin/time --output=lengths/timealt12_p010_27.txt --format="\t%E real,\t%U user,\t%S sys" ./sl12_cl_alt_noselfcheck.exe 200001 0.50 1.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/timealt12_p030_27.txt --format="\t%E real,\t%U user,\t%S sys" ./sl12_cl_alt_noselfcheck.exe 200001 0.50 3.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/timealt12_p050_27.txt --format="\t%E real,\t%U user,\t%S sys" ./sl12_cl_alt_noselfcheck.exe 200001 0.50 5.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/timealt12_p070_27.txt --format="\t%E real,\t%U user,\t%S sys" ./sl12_cl_alt_noselfcheck.exe 200001 0.50 7.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/timealt12_p080_27.txt --format="\t%E real,\t%U user,\t%S sys" ./sl12_cl_alt_noselfcheck.exe 200001 0.50 8.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/timealt12_p100_27.txt --format="\t%E real,\t%U user,\t%S sys" ./sl12_cl_alt_noselfcheck.exe 200001 0.50 10.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/timealt12_p120_27.txt --format="\t%E real,\t%U user,\t%S sys" ./sl12_cl_alt_noselfcheck.exe 200001 0.50 12.0 1.26610367 > /dev/null &
/usr/bin/time --output=lengths/timealt12_p150_27.txt --format="\t%E real,\t%U user,\t%S sys" ./sl12_cl_alt_noselfcheck.exe 200001 0.50 15.0 1.26610367 > /dev/null &
