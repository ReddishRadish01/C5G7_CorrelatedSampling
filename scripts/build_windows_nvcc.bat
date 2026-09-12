@echo off
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0build_windows_nvcc.ps1" %*
