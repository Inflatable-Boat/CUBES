#!/bin/sh
/usr/bin/time -o times/vtemp_20_a1.16 ./vtemp c 20 300001 0.50 10 1.57079632679489 > logs/vtemp_20_a16 &
/usr/bin/time -o times/vtemp_20_a1.27 ./vtemp c 20 300001 0.50 10 1.47062890563334 > logs/vtemp_20_a27 &
/usr/bin/time -o times/vtemp_20_a1.37 ./vtemp c 20 300001 0.50 10 1.36943840600457 > logs/vtemp_20_a37 &
/usr/bin/time -o times/vtemp_20_a1.47 ./vtemp c 20 300001 0.50 10 1.26610367277950 > logs/vtemp_20_a47 &
/usr/bin/time -o times/vtemp_20_a1.57 ./vtemp c 20 300001 0.50 10 1.15927948072741 > logs/vtemp_20_a57 &
