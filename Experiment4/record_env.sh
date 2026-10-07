#!/usr/bin/env bash
# Records the facts needed for the manuscript caption. Works in MSYS2 (Windows) and Linux.
# Usage: ./record_env.sh > results/env.txt
echo "--- date ---"; date
echo "--- CPU ---"
if command -v lscpu >/dev/null 2>&1; then
  lscpu
elif command -v powershell.exe >/dev/null 2>&1; then
  powershell.exe -NoProfile -Command "Get-CimInstance Win32_Processor | Format-List Name,NumberOfCores,NumberOfLogicalProcessors,MaxClockSpeed"
fi
echo "--- RAM ---"
if command -v free >/dev/null 2>&1 && [ -r /proc/meminfo ]; then
  grep MemTotal /proc/meminfo
elif command -v powershell.exe >/dev/null 2>&1; then
  powershell.exe -NoProfile -Command "[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB,1)"
fi
echo "--- OS ---"
uname -a
if command -v powershell.exe >/dev/null 2>&1; then
  powershell.exe -NoProfile -Command "(Get-CimInstance Win32_OperatingSystem).Caption"
  echo "--- Windows power plan ---"
  powershell.exe -NoProfile -Command "powercfg /getactivescheme"
fi
echo "--- compiler ---"; g++ --version | head -1
echo "--- OpenMP macro ---"; echo | g++ -fopenmp -dM -E - | grep -i _OPENMP
echo "--- GMP package ---"
if command -v pacman >/dev/null 2>&1; then pacman -Q | grep -i gmp
elif command -v dpkg >/dev/null 2>&1; then dpkg -l | grep -i gmp; fi
echo "--- OMP env ---"
echo "OMP_PROC_BIND=${OMP_PROC_BIND:-<unset>}  OMP_PLACES=${OMP_PLACES:-<unset>}  OMP_NUM_THREADS=${OMP_NUM_THREADS:-<unset>}"
