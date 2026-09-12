@echo off
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_critical_rod_sweep_windows_nvcc.ps1" %*
