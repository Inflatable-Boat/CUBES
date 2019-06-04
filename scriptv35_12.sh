#!/bin/sh
/usr/bin/time -o times/v35_12_a57 ./v35_sim.exe c 12 0.3 500000 2000 1.57079632679489 > logs/v35_12_a57 &
/usr/bin/time -o times/v35_12_a47 ./v35_sim.exe c 12 0.3 500000 2000 1.47062890563334 > logs/v35_12_a47 &
/usr/bin/time -o times/v35_12_a37 ./v35_sim.exe c 12 0.3 500000 2000 1.36943840600457 > logs/v35_12_a37 &
/usr/bin/time -o times/v35_12_a27 ./v35_sim.exe c 12 0.3 500000 2000 1.26610367277950 > logs/v35_12_a27 &
/usr/bin/time -o times/v35_12_a16 ./v35_sim.exe c 12 0.3 500000 2000 1.15927948072741 > logs/v35_12_a16 &
