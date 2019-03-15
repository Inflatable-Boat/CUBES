#!/bin/sh

time -o times/vumbr_20p5.0 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.57.poly vumbr_20p5.0 200000 200 5.0 1.57079632679 > logs/vumbr_20p5.0 &
time -o times/vumbr_20p5.2 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.57.poly vumbr_20p5.2 200000 200 5.2 1.57079632679 > logs/vumbr_20p5.2 &
time -o times/vumbr_20p5.4 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.57.poly vumbr_20p5.4 200000 200 5.4 1.57079632679 > logs/vumbr_20p5.4 &
time -o times/vumbr_20p5.6 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.57.poly vumbr_20p5.6 200000 200 5.6 1.57079632679 > logs/vumbr_20p5.6 &
time -o times/vumbr_20p5.8 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.57.poly vumbr_20p5.8 200000 200 5.8 1.57079632679 > logs/vumbr_20p5.8 &
time -o times/vumbr_20p6.0 -f "%E\t%U\t%S,\t%P\n" ./v21_sim.exe r v19_20_a1.57.poly vumbr_20p6.0 200000 200 6.0 1.57079632679 > logs/vumbr_20p6.0 &