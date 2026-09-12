@echo off
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_critical_rod_newton_cs_windows_nvcc.ps1" %*
exit /b %ERRORLEVEL%
