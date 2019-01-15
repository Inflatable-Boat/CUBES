#!/bin/sh
/usr/bin/time -o times/v3_16_a57 -f "%E, %P\n" ./v3_read c 16 2000001 0.50 1.57079632679489 > logs/v3_16_a57 &
/usr/bin/time -o times/v3_16_a47 -f "%E, %P\n" ./v3_read c 16 2000001 0.50 1.47062890563334 > logs/v3_16_a47 &
/usr/bin/time -o times/v3_16_a37 -f "%E, %P\n" ./v3_read c 16 2000001 0.50 1.36943840600457 > logs/v3_16_a37 &
/usr/bin/time -o times/v3_16_a27 -f "%E, %P\n" ./v3_read c 16 2000001 0.50 1.26610367277950 > logs/v3_16_a27 &
/usr/bin/time -o times/v3_16_a16 -f "%E, %P\n" ./v3_read c 16 2000001 0.50 1.15927948072741 > logs/v3_16_a16 &
