#!/bin/sh
/usr/bin/time -o times/v13_20_a57 -f "%E, %P\n" ./v13_crush_slowly c 20 1000001 1.57079632679489 > logs/v13_20_a57 &
/usr/bin/time -o times/v13_20_a47 -f "%E, %P\n" ./v13_crush_slowly c 20 1000001 1.47062890563334 > logs/v13_20_a47 &
/usr/bin/time -o times/v13_20_a37 -f "%E, %P\n" ./v13_crush_slowly c 20 1000001 1.36943840600457 > logs/v13_20_a37 &
/usr/bin/time -o times/v13_20_a27 -f "%E, %P\n" ./v13_crush_slowly c 20 1000001 1.26610367277950 > logs/v13_20_a27 &
/usr/bin/time -o times/v13_20_a16 -f "%E, %P\n" ./v13_crush_slowly c 20 1000001 1.15927948072741 > logs/v13_20_a16 &
