@echo off
REM UniPPI-OWL v1.8 Local Loopback Launcher (Windows)
echo ======================================================================
echo  UniPPI-OWL v1.8 Local App: Opening http://127.0.0.1:8766/
echo  Data Privacy Lock: 100%% Local Execution (Close window to stop)
echo ======================================================================
start http://127.0.0.1:8766/
python -m http.server 8766 --bind 127.0.0.1
