@echo off
for /f "usebackq tokens=1-6 delims=," %%a in ("parameters.csv") do (
      Set "_var1=%%a"
      Set "_var2=%%b"
      Set "_var3=%%c"
      Set "_var4=%%d"
      Set "_var5=%%e"
      Set "_var6=%%f"
)
python PTM_predictor.py %_var1% %_var2% %_var3% %_var4% %_var5% %_var6%
pause
