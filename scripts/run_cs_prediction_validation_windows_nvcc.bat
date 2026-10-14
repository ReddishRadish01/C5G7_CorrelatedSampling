@echo off
where pwsh.exe >nul 2>&1
if errorlevel 1 goto windows_powershell
pwsh.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_cs_prediction_validation_windows_nvcc.ps1" %*
exit /b %ERRORLEVEL%

:windows_powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_cs_prediction_validation_windows_nvcc.ps1" %*
exit /b %ERRORLEVEL%
