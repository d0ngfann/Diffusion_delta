@echo off
cd /d %~dp0

echo === Simulation p2q4 ===
python sweep_simulation/sweep_p2_vs_thrshld.py --p 2 --q 4 --n-pairs 500 --n-workers 14
echo === Simulation p3q2 ===
python sweep_simulation/sweep_p2_vs_thrshld.py --p 3 --q 2 --n-pairs 500 --n-workers 14

echo === ALL DONE ===
pause
