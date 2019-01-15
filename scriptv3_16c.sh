#!/bin/sh
/usr/bin/time -o times/v3_16_a57c -f "%E, %P\n" ./order_makerv3 16 2000001 0.50 1.0 1.57079632679489 > logs/v3_16_a57c &
/usr/bin/time -o times/v3_16_a47c -f "%E, %P\n" ./order_makerv3 16 2000001 0.50 1.0 1.47062890563334 > logs/v3_16_a47c &
/usr/bin/time -o times/v3_16_a37c -f "%E, %P\n" ./order_makerv3 16 2000001 0.50 1.0 1.36943840600457 > logs/v3_16_a37c &
/usr/bin/time -o times/v3_16_a27c -f "%E, %P\n" ./order_makerv3 16 2000001 0.50 1.0 1.26610367277950 > logs/v3_16_a27c &
/usr/bin/time -o times/v3_16_a16c -f "%E, %P\n" ./order_makerv3 16 2000001 0.50 1.0 1.15927948072741 > logs/v3_16_a16c &
