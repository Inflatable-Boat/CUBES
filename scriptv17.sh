#!/bin/sh
time -o times/v17_12_pf0.40_p99.0_a1.27_os_t0010_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 99.0 1.26610367277950 s 10 npt 0.005 > logs/v17_12_pf0.40_p99.0_a1.27_os_t0010_npt_c0.005 &
time -o times/v17_12_pf0.40_p12.0_a1.27_os_t0010_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 12.0 1.26610367277950 s 10 npt 0.005 > logs/v17_12_pf0.40_p12.0_a1.27_os_t0010_npt_c0.005 &
time -o times/v17_12_pf0.40_p99.0_a1.27_os_t0100_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 99.0 1.26610367277950 s 100 npt 0.005 > logs/v17_12_pf0.40_p99.0_a1.27_os_t0100_npt_c0.005 &
time -o times/v17_12_pf0.40_p12.0_a1.27_os_t0100_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 12.0 1.26610367277950 s 100 npt 0.005 > logs/v17_12_pf0.40_p12.0_a1.27_os_t0100_npt_c0.005 &
time -o times/v17_12_pf0.40_p99.0_a1.27_os_t0010_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 99.0 1.26610367277950 s 10 npt 0.001 > logs/v17_12_pf0.40_p99.0_a1.27_os_t0010_npt_c0.001 &
time -o times/v17_12_pf0.40_p12.0_a1.27_os_t0010_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 12.0 1.26610367277950 s 10 npt 0.001 > logs/v17_12_pf0.40_p12.0_a1.27_os_t0010_npt_c0.001 &
time -o times/v17_12_pf0.40_p99.0_a1.27_os_t0100_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 99.0 1.26610367277950 s 100 npt 0.001 > logs/v17_12_pf0.40_p99.0_a1.27_os_t0100_npt_c0.001 &
time -o times/v17_12_pf0.40_p12.0_a1.27_os_t0100_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 12.0 1.26610367277950 s 100 npt 0.001 > logs/v17_12_pf0.40_p12.0_a1.27_os_t0100_npt_c0.001 &

time -o times/v17_12_pf0.40_p99.0_a1.27_ot_t0010_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 99.0 1.26610367277950 t 10 npt 0.005 > logs/v17_12_pf0.40_p99.0_a1.27_ot_t0010_npt_c0.005 &
time -o times/v17_12_pf0.40_p12.0_a1.27_ot_t0010_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 12.0 1.26610367277950 t 10 npt 0.005 > logs/v17_12_pf0.40_p12.0_a1.27_ot_t0010_npt_c0.005 &
time -o times/v17_12_pf0.40_p99.0_a1.27_ot_t0100_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 99.0 1.26610367277950 t 100 npt 0.005 > logs/v17_12_pf0.40_p99.0_a1.27_ot_t0100_npt_c0.005 &
time -o times/v17_12_pf0.40_p12.0_a1.27_ot_t0100_npt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 12.0 1.26610367277950 t 100 npt 0.005 > logs/v17_12_pf0.40_p12.0_a1.27_ot_t0100_npt_c0.005 &
time -o times/v17_12_pf0.40_p99.0_a1.27_ot_t0010_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 99.0 1.26610367277950 t 10 npt 0.001 > logs/v17_12_pf0.40_p99.0_a1.27_ot_t0010_npt_c0.001 &
time -o times/v17_12_pf0.40_p12.0_a1.27_ot_t0010_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 12.0 1.26610367277950 t 10 npt 0.001 > logs/v17_12_pf0.40_p12.0_a1.27_ot_t0010_npt_c0.001 &
time -o times/v17_12_pf0.40_p99.0_a1.27_ot_t0100_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 99.0 1.26610367277950 t 100 npt 0.001 > logs/v17_12_pf0.40_p99.0_a1.27_ot_t0100_npt_c0.001 &
time -o times/v17_12_pf0.40_p12.0_a1.27_ot_t0100_npt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.40 500000 100 12.0 1.26610367277950 t 100 npt 0.001 > logs/v17_12_pf0.40_p12.0_a1.27_ot_t0100_npt_c0.001 &

time -o times/v17_12_pf0.50_p99.0_a1.27_os_t0010_nvt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.50 500000 100 99.0 1.26610367277950 s 10 nvt 0.005 > logs/v17_12_pf0.50_p99.0_a1.27_os_t0010_nvt_c0.005 &
time -o times/v17_12_pf0.50_p99.0_a1.27_os_t0010_nvt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.50 500000 100 99.0 1.26610367277950 s 10 nvt 0.001 > logs/v17_12_pf0.50_p99.0_a1.27_os_t0010_nvt_c0.001 &
time -o times/v17_12_pf0.50_p99.0_a1.27_os_t0100_nvt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.50 500000 100 99.0 1.26610367277950 s 100 nvt 0.005 > logs/v17_12_pf0.50_p99.0_a1.27_os_t0100_nvt_c0.005 &
time -o times/v17_12_pf0.50_p99.0_a1.27_os_t0100_nvt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.50 500000 100 99.0 1.26610367277950 s 100 nvt 0.001 > logs/v17_12_pf0.50_p99.0_a1.27_os_t0100_nvt_c0.001 &

time -o times/v17_12_pf0.50_p99.0_a1.27_ot_t0010_nvt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.50 500000 100 99.0 1.26610367277950 t 10 nvt 0.005 > logs/v17_12_pf0.50_p99.0_a1.27_ot_t0010_nvt_c0.005 &
time -o times/v17_12_pf0.50_p99.0_a1.27_ot_t0010_nvt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.50 500000 100 99.0 1.26610367277950 t 10 nvt 0.001 > logs/v17_12_pf0.50_p99.0_a1.27_ot_t0010_nvt_c0.001 &
time -o times/v17_12_pf0.50_p99.0_a1.27_ot_t0100_nvt_c0.005 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.50 500000 100 99.0 1.26610367277950 t 100 nvt 0.005 > logs/v17_12_pf0.50_p99.0_a1.27_ot_t0100_nvt_c0.005 &
time -o times/v17_12_pf0.50_p99.0_a1.27_ot_t0100_nvt_c0.001 -f "%E, %P" ./v17_sim_umbrella.exe c 12 0.50 500000 100 99.0 1.26610367277950 t 100 nvt 0.001 > logs/v17_12_pf0.50_p99.0_a1.27_ot_t0100_nvt_c0.001 &
