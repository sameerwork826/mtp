@echo off
:: Script to run LAMMPS simulations with different seeds on Windows
:: Place this script in the same directory as your LAMMPS input file and data file

:: Original LAMMPS input file
set ORIGINAL_FILE=in_poly.lammps

:: Create and modify the first input file with seed 12345
set NEW_FILE=in.poly.seed1.lammps
copy %ORIGINAL_FILE% %NEW_FILE% > nul
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'fix\s+1\s+all\s+langevin\s+\${temp}\s+\${temp}\s+1\.0\s+4', 'fix             1 all langevin ${temp} ${temp} 1.0 12345' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'compute\s+msd_12_10_4', 'compute         msd_12_10_seed1' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'print \"Time msd_12_10_4\"', 'print \"Time msd_12_10_seed1\"' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'file msd_12_10_4.dat', 'file msd_12_10_seed1.dat' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'c_msd_12_10_4\[4\]', 'c_msd_12_10_seed1[4]' } | Set-Content %NEW_FILE%"
echo Created input file %NEW_FILE% with seed 12345

:: Create and modify the second input file with seed 67890
set NEW_FILE=in.poly.seed2.lammps
copy %ORIGINAL_FILE% %NEW_FILE% > nul
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'fix\s+1\s+all\s+langevin\s+\${temp}\s+\${temp}\s+1\.0\s+4', 'fix             1 all langevin ${temp} ${temp} 1.0 67890' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'compute\s+msd_12_10_4', 'compute         msd_12_10_seed2' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'print \"Time msd_12_10_4\"', 'print \"Time msd_12_10_seed2\"' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'file msd_12_10_4.dat', 'file msd_12_10_seed2.dat' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'c_msd_12_10_4\[4\]', 'c_msd_12_10_seed2[4]' } | Set-Content %NEW_FILE%"
echo Created input file %NEW_FILE% with seed 67890

:: Create and modify the third input file with seed 24680
set NEW_FILE=in.poly.seed3.lammps
copy %ORIGINAL_FILE% %NEW_FILE% > nul
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'fix\s+1\s+all\s+langevin\s+\${temp}\s+\${temp}\s+1\.0\s+4', 'fix             1 all langevin ${temp} ${temp} 1.0 24680' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'compute\s+msd_12_10_4', 'compute         msd_12_10_seed3' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'print \"Time msd_12_10_4\"', 'print \"Time msd_12_10_seed3\"' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'file msd_12_10_4.dat', 'file msd_12_10_seed3.dat' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'c_msd_12_10_4\[4\]', 'c_msd_12_10_seed3[4]' } | Set-Content %NEW_FILE%"
echo Created input file %NEW_FILE% with seed 24680

:: Create and modify the fourth input file with seed 13579
set NEW_FILE=in.poly.seed4.lammps
copy %ORIGINAL_FILE% %NEW_FILE% > nul
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'fix\s+1\s+all\s+langevin\s+\${temp}\s+\${temp}\s+1\.0\s+4', 'fix             1 all langevin ${temp} ${temp} 1.0 13579' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'compute\s+msd_12_10_4', 'compute         msd_12_10_seed4' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'print \"Time msd_12_10_4\"', 'print \"Time msd_12_10_seed4\"' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'file msd_12_10_4.dat', 'file msd_12_10_seed4.dat' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'c_msd_12_10_4\[4\]', 'c_msd_12_10_seed4[4]' } | Set-Content %NEW_FILE%"
echo Created input file %NEW_FILE% with seed 13579

:: Launch each simulation in a separate CMD window
echo Starting simulations...
start cmd /k "lmp -in in.poly.seed1.lammps && echo Simulation 1 completed"
start cmd /k "lmp -in in.poly.seed2.lammps && echo Simulation 2 completed"
start cmd /k "lmp -in in.poly.seed3.lammps && echo Simulation 3 completed"
start cmd /k "lmp -in in.poly.seed4.lammps && echo Simulation 4 completed"

echo All simulations launched!