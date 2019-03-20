#!/bin/sh
time -o times/v22_vumbr_20p5.0 -f "%E\t%U\t%S,\t%P\n" ./v22_coords_to_order.exe vumbr_20p5.0 200 1000 b s > logs/v22_vumbr_20p5.0 &
time -o times/v22_vumbr_20p5.2 -f "%E\t%U\t%S,\t%P\n" ./v22_coords_to_order.exe vumbr_20p5.2 200 1000 b s > logs/v22_vumbr_20p5.2 &
time -o times/v22_vumbr_20p5.4 -f "%E\t%U\t%S,\t%P\n" ./v22_coords_to_order.exe vumbr_20p5.4 200 1000 b s > logs/v22_vumbr_20p5.4 &
time -o times/v22_vumbr_20p5.6 -f "%E\t%U\t%S,\t%P\n" ./v22_coords_to_order.exe vumbr_20p5.6 200 1000 b s > logs/v22_vumbr_20p5.6 &
time -o times/v22_vumbr_20p5.8 -f "%E\t%U\t%S,\t%P\n" ./v22_coords_to_order.exe vumbr_20p5.8 200 1000 b s > logs/v22_vumbr_20p5.8 &
time -o times/v22_vumbr_20p6.0 -f "%E\t%U\t%S,\t%P\n" ./v22_coords_to_order.exe vumbr_20p6.0 200 1000 b s > logs/v22_vumbr_20p6.0 &