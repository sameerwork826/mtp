@echo off
:: Script to run LAMMPS simulations with different seeds on Windows
:: Place this script in the same directory as your LAMMPS input file and data file
:: Original LAMMPS input file
set ORIGINAL_FILE=C:\Users\nande\Desktop\mtp\in.lammps  

:: Create and modify the first input file with seed 12345
set NEW_FILE=in.poly_12.seed1.lammps
copy %ORIGINAL_FILE% %NEW_FILE% > nul
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'fix\s+1\s+all\s+langevin\s+\${temp}\s+\${temp}\s+1\.0\s+4', 'fix 1 all langevin ${temp} ${temp} 1.0 12345' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'file msd_all.dat', 'file msd_all_12_1.dat' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'fix\s+rg_avg_all\s+all\s+ave/time\s+100\s+1\s+100\s+c_rg_all\s+file\s+rg_all.dat', 'fix rg_avg_all all ave/time 100 1 100 c_rg_all file rg_all_12_1.dat' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'dump\s+prod_dump\s+all\s+custom\s+1000\s+production.lammpstrj', 'dump prod_dump all custom 1000 production_12_1.lammpstrj' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'dump\s+eq_dump\s+all\s+custom\s+1000\s+equilibration.lammpstrj', 'dump eq_dump all custom 1000 equilibration_12_1.lammpstrj' } | Set-Content %NEW_FILE%"
echo Created input file %NEW_FILE% with seed 12345            

:: Create and modify the second input file with seed 67890
set NEW_FILE=in.poly_12.seed2.lammps
copy %ORIGINAL_FILE% %NEW_FILE% > nul
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'fix\s+1\s+all\s+langevin\s+\${temp}\s+\${temp}\s+1\.0\s+4', 'fix 1 all langevin ${temp} ${temp} 1.0 67890' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'file msd_all.dat', 'file msd_all_12_2.dat' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'fix\s+rg_avg_all\s+all\s+ave/time\s+100\s+1\s+100\s+c_rg_all\s+file\s+rg_all.dat', 'fix rg_avg_all all ave/time 100 1 100 c_rg_all file rg_all_12_2.dat' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'dump\s+prod_dump\s+all\s+custom\s+1000\s+production.lammpstrj', 'dump prod_dump all custom 1000 production_12_2.lammpstrj' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'dump\s+eq_dump\s+all\s+custom\s+1000\s+equilibration.lammpstrj', 'dump eq_dump all custom 1000 equilibration_12_2.lammpstrj' } | Set-Content %NEW_FILE%"
echo Created input file %NEW_FILE% with seed 67890

:: Create and modify the third input file with seed 12680
set NEW_FILE=in.poly_12.seed3.lammps
copy %ORIGINAL_FILE% %NEW_FILE% > nul
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'fix\s+1\s+all\s+langevin\s+\${temp}\s+\${temp}\s+1\.0\s+4', 'fix 1 all langevin ${temp} ${temp} 1.0 12680' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'file msd_all.dat', 'file msd_all_12_3.dat' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'fix\s+rg_avg_all\s+all\s+ave/time\s+100\s+1\s+100\s+c_rg_all\s+file\s+rg_all.dat', 'fix rg_avg_all all ave/time 100 1 100 c_rg_all file rg_all_12_3.dat' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'dump\s+prod_dump\s+all\s+custom\s+1000\s+production.lammpstrj', 'dump prod_dump all custom 1000 production_12_3.lammpstrj' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'dump\s+eq_dump\s+all\s+custom\s+1000\s+equilibration.lammpstrj', 'dump eq_dump all custom 1000 equilibration_12_3.lammpstrj' } | Set-Content %NEW_FILE%"
echo Created input file %NEW_FILE% with seed 12680

:: Create and modify the fourth input file with seed 13579
set NEW_FILE=in.poly_12.seed4.lammps
copy %ORIGINAL_FILE% %NEW_FILE% > nul
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'fix\s+1\s+all\s+langevin\s+\${temp}\s+\${temp}\s+1\.0\s+4', 'fix 1 all langevin ${temp} ${temp} 1.0 13579' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'file msd_all.dat', 'file msd_all_12_4.dat' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'fix\s+rg_avg_all\s+all\s+ave/time\s+100\s+1\s+100\s+c_rg_all\s+file\s+rg_all.dat', 'fix rg_avg_all all ave/time 100 1 100 c_rg_all file rg_all_12_4.dat' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'dump\s+prod_dump\s+all\s+custom\s+1000\s+production.lammpstrj', 'dump prod_dump all custom 1000 production_12_4.lammpstrj' } | Set-Content %NEW_FILE%"
powershell -Command "(Get-Content %NEW_FILE%) | ForEach-Object { $_ -replace 'dump\s+eq_dump\s+all\s+custom\s+1000\s+equilibration.lammpstrj', 'dump eq_dump all custom 1000 equilibration_12_4.lammpstrj' } | Set-Content %NEW_FILE%"
echo Created input file %NEW_FILE% with seed 13579

:: Launch each simulation in a separate CMD window
echo Starting simulations...
start cmd /k "lmp -in in.poly_12.seed1.lammps && echo Simulation 1 completed"
start cmd /k "lmp -in in.poly_12.seed2.lammps && echo Simulation 2 completed"
start cmd /k "lmp -in in.poly_12.seed3.lammps && echo Simulation 3 completed"
start cmd /k "lmp -in in.poly_12.seed4.lammps && echo Simulation 4 completed"
echo All simulations launched!