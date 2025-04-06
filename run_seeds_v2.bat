@echo off
setlocal EnableDelayedExpansion

:: Script to run LAMMPS simulations with different data files and seeds on Windows
:: Place this script in the same directory as your LAMMPS executable

:: Original LAMMPS input file
set ORIGINAL_FILE=C:\Users\nande\Desktop\mtp\in_poly.lammps

:: Define the base path for data files
set DATA_PATH=C:\Users\nande\Desktop\mtp\case_1_only_dn\data_files

:: Define random seeds to use
set SEED1=12345
set SEED2=67890
set SEED3=12680
set SEED4=13579

:: Array of polymer sizes to process
set SIZES=12 16 20 24 28 32

:: Process each polymer size
for %%s in (%SIZES%) do (
    echo Processing polymer size %%s...
    
    :: Process each seed for this polymer size
    for /L %%i in (1,1,4) do (
        :: Set the seed based on the current iteration
        if %%i==1 set CURRENT_SEED=%SEED1%
        if %%i==2 set CURRENT_SEED=%SEED2%
        if %%i==3 set CURRENT_SEED=%SEED3%
        if %%i==4 set CURRENT_SEED=%SEED4%
        
        :: Create the new input file name
        set NEW_FILE=in.poly_%%s.seed%%i.lammps
        
        :: Copy the original file to create a new input file
        copy "%ORIGINAL_FILE%" "!NEW_FILE!" > nul
        
        :: Update the data file path - using delayed expansion for variables
        powershell -Command "(Get-Content '!NEW_FILE!') -replace 'read_data\s+.*\\polymer_system_\d+\.data', 'read_data %DATA_PATH%\polymer_system_%%s.data' | Set-Content '!NEW_FILE!'"
        
        :: Update the langevin thermostat fix with the new seed
        powershell -Command "(Get-Content '!NEW_FILE!') -replace 'fix\s+1\s+all\s+langevin\s+\${temp}\s+\${temp}\s+\d+\.\d+\s+\d+', 'fix 1 all langevin ${temp} ${temp} 2.0 !CURRENT_SEED!' | Set-Content '!NEW_FILE!'"
        
        :: Update output filenames with polymer size and seed number
        powershell -Command "(Get-Content '!NEW_FILE!') -replace 'file\s+msd_all\.dat', 'file msd_all_%%s_%%i.dat' | Set-Content '!NEW_FILE!'"
        powershell -Command "(Get-Content '!NEW_FILE!') -replace 'fix\s+rg_avg_all\s+all\s+ave/time\s+\d+\s+\d+\s+\d+\s+c_rg_all\s+file\s+rg_all\.dat', 'fix rg_avg_all all ave/time 100 1 100 c_rg_all file rg_all_%%s_%%i.dat' | Set-Content '!NEW_FILE!'"
        powershell -Command "(Get-Content '!NEW_FILE!') -replace 'dump\s+prod_dump\s+all\s+custom\s+\d+\s+production\.lammpstrj', 'dump prod_dump all custom 1000 production_%%s_%%i.lammpstrj' | Set-Content '!NEW_FILE!'"
        powershell -Command "(Get-Content '!NEW_FILE!') -replace 'dump\s+eq_dump\s+all\s+custom\s+\d+\s+equilibration\.lammpstrj', 'dump eq_dump all custom 1000 equilibration_%%s_%%i.lammpstrj' | Set-Content '!NEW_FILE!'"
        
        echo Created input file !NEW_FILE! with seed !CURRENT_SEED! for polymer size %%s
    )
)

echo All input files created!

:: Ask whether to launch simulations
set /p LAUNCH_SIMS=Do you want to launch all simulations now? (Y/N): 

if /i "%LAUNCH_SIMS%"=="Y" (
    echo Starting simulations...
    
    :: Launch each simulation in a separate CMD window
    for %%s in (%SIZES%) do (
        for /L %%i in (1,1,4) do (
            start cmd /k "lmp -in in.poly_%%s.seed%%i.lammps && echo Simulation for polymer size %%s with seed %%i completed"
            :: Add a short delay to prevent overwhelming the system
            timeout /t 2 /nobreak > nul
        )
    )
    
    echo All simulations launched!
) else (
    echo Simulation files are ready to run manually.
)