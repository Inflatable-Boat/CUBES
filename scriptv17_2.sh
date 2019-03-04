#!/bin/sh
time -o times/v17_12_pf0.45_p08.0_a1.27_os_t0050_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 8.0 1.26610367277950 s 50  npt 0.005 > logs/v17_12_pf0.45_p08.0_a1.27_os_t0050_npt_c0.005 &
time -o times/v17_12_pf0.45_p07.0_a1.27_os_t0050_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 7.0 1.26610367277950 s 50  npt 0.005 > logs/v17_12_pf0.45_p07.0_a1.27_os_t0050_npt_c0.005 &
time -o times/v17_12_pf0.45_p08.0_a1.27_os_t0200_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 8.0 1.26610367277950 s 200 npt 0.005 > logs/v17_12_pf0.45_p08.0_a1.27_os_t0200_npt_c0.005 &
time -o times/v17_12_pf0.45_p07.0_a1.27_os_t0200_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 7.0 1.26610367277950 s 200 npt 0.005 > logs/v17_12_pf0.45_p07.0_a1.27_os_t0200_npt_c0.005 &

time -o times/v17_12_pf0.45_p08.0_a1.27_ot_t0050_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 8.0 1.26610367277950 t 50  npt 0.005 > logs/v17_12_pf0.45_p08.0_a1.27_ot_t0050_npt_c0.005 &
time -o times/v17_12_pf0.45_p07.0_a1.27_ot_t0050_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 7.0 1.26610367277950 t 50  npt 0.005 > logs/v17_12_pf0.45_p07.0_a1.27_ot_t0050_npt_c0.005 &
time -o times/v17_12_pf0.45_p08.0_a1.27_ot_t0200_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 8.0 1.26610367277950 t 200 npt 0.005 > logs/v17_12_pf0.45_p08.0_a1.27_ot_t0200_npt_c0.005 &
time -o times/v17_12_pf0.45_p07.0_a1.27_ot_t0200_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 7.0 1.26610367277950 t 200 npt 0.005 > logs/v17_12_pf0.45_p07.0_a1.27_ot_t0200_npt_c0.005 &

time -o times/v17_12_pf0.45_p08.0_a1.27_os_t0050_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 8.0 1.26610367277950 s 50  npt 0.001 > logs/v17_12_pf0.45_p08.0_a1.27_os_t0050_npt_c0.001 &
time -o times/v17_12_pf0.45_p07.0_a1.27_os_t0050_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 7.0 1.26610367277950 s 50  npt 0.001 > logs/v17_12_pf0.45_p07.0_a1.27_os_t0050_npt_c0.001 &
time -o times/v17_12_pf0.45_p08.0_a1.27_os_t0200_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 8.0 1.26610367277950 s 200 npt 0.001 > logs/v17_12_pf0.45_p08.0_a1.27_os_t0200_npt_c0.001 &
time -o times/v17_12_pf0.45_p07.0_a1.27_os_t0200_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 7.0 1.26610367277950 s 200 npt 0.001 > logs/v17_12_pf0.45_p07.0_a1.27_os_t0200_npt_c0.001 &

time -o times/v17_12_pf0.45_p08.0_a1.27_ot_t0050_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 8.0 1.26610367277950 t 50  npt 0.001 > logs/v17_12_pf0.45_p08.0_a1.27_ot_t0050_npt_c0.001 &
time -o times/v17_12_pf0.45_p07.0_a1.27_ot_t0050_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 7.0 1.26610367277950 t 50  npt 0.001 > logs/v17_12_pf0.45_p07.0_a1.27_ot_t0050_npt_c0.001 &
time -o times/v17_12_pf0.45_p08.0_a1.27_ot_t0200_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 8.0 1.26610367277950 t 200 npt 0.001 > logs/v17_12_pf0.45_p08.0_a1.27_ot_t0200_npt_c0.001 &
time -o times/v17_12_pf0.45_p07.0_a1.27_ot_t0200_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.45 500000 100 7.0 1.26610367277950 t 200 npt 0.001 > logs/v17_12_pf0.45_p07.0_a1.27_ot_t0200_npt_c0.001 &
