#!/bin/sh
/usr/bin/time -o times/v5_10_p10a57 -f "%E, %P\n" ./v5_read_g_of_r c 10 200001 0.50 10.0 1.57079632679489 > logs/v5_10_p10a57 &
/usr/bin/time -o times/v5_10_p10a47 -f "%E, %P\n" ./v5_read_g_of_r c 10 200001 0.50 10.0 1.47062890563334 > logs/v5_10_p10a47 &
/usr/bin/time -o times/v5_10_p10a37 -f "%E, %P\n" ./v5_read_g_of_r c 10 200001 0.50 10.0 1.36943840600457 > logs/v5_10_p10a37 &
/usr/bin/time -o times/v5_10_p10a27 -f "%E, %P\n" ./v5_read_g_of_r c 10 200001 0.50 10.0 1.26610367277950 > logs/v5_10_p10a27 &
/usr/bin/time -o times/v5_10_p10a16 -f "%E, %P\n" ./v5_read_g_of_r c 10 200001 0.50 10.0 1.15927948072741 > logs/v5_10_p10a16 &

/usr/bin/time -o times/v5_10_p05a57 -f "%E, %P\n" ./v5_read_g_of_r c 10 200001 0.50 5.0 1.57079632679489 > logs/v5_10_p05a57 &
/usr/bin/time -o times/v5_10_p05a47 -f "%E, %P\n" ./v5_read_g_of_r c 10 200001 0.50 5.0 1.47062890563334 > logs/v5_10_p05a47 &
/usr/bin/time -o times/v5_10_p05a37 -f "%E, %P\n" ./v5_read_g_of_r c 10 200001 0.50 5.0 1.36943840600457 > logs/v5_10_p05a37 &
/usr/bin/time -o times/v5_10_p05a27 -f "%E, %P\n" ./v5_read_g_of_r c 10 200001 0.50 5.0 1.26610367277950 > logs/v5_10_p05a27 &
/usr/bin/time -o times/v5_10_p05a16 -f "%E, %P\n" ./v5_read_g_of_r c 10 200001 0.50 5.0 1.15927948072741 > logs/v5_10_p05a16 &
