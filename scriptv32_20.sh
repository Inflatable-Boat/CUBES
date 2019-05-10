#!/bin/sh
/usr/bin/time -o times/v32_20_a57 ./v32_sim_crush.exe c 20 0.1 2000001 1000 1.57079632679489 > logs/v32_20_a57 &
/usr/bin/time -o times/v32_20_a47 ./v32_sim_crush.exe c 20 0.1 2000001 1000 1.47062890563334 > logs/v32_20_a47 &
/usr/bin/time -o times/v32_20_a37 ./v32_sim_crush.exe c 20 0.1 2000001 1000 1.36943840600457 > logs/v32_20_a37 &
/usr/bin/time -o times/v32_20_a27 ./v32_sim_crush.exe c 20 0.1 2000001 1000 1.26610367277950 > logs/v32_20_a27 &
/usr/bin/time -o times/v32_20_a16 ./v32_sim_crush.exe c 20 0.1 2000001 1000 1.15927948072741 > logs/v32_20_a16 &
