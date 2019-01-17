#!/bin/sh
/usr/bin/time -o times/v4_20_a57 -f "%E, %P\n" ./v4_read c 20 2000001 1.57079632679489 > logs/v4_20_a57 &
/usr/bin/time -o times/v4_20_a47 -f "%E, %P\n" ./v4_read c 20 2000001 1.47062890563334 > logs/v4_20_a47 &
/usr/bin/time -o times/v4_20_a37 -f "%E, %P\n" ./v4_read c 20 2000001 1.36943840600457 > logs/v4_20_a37 &
/usr/bin/time -o times/v4_20_a27 -f "%E, %P\n" ./v4_read c 20 2000001 1.26610367277950 > logs/v4_20_a27 &
/usr/bin/time -o times/v4_20_a16 -f "%E, %P\n" ./v4_read c 20 2000001 1.15927948072741 > logs/v4_20_a16 &
