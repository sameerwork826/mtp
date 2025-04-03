@echo off
setlocal enabledelayedexpansion

:: Configuration
set "baseDir=C:\Users\nande\Desktop\mtp\case_3_og"
set "dataFile=%baseDir%\data_files\ps_c_5_12.data"
set "resultsDir=%baseDir%\results\12"
set "lammpsExec=lmp.exe"  :: Change to your LAMMPS executable

:: Seeds to run
set seeds=12345 54321 98765 13579

:: Create results directory
if not exist "%resultsDir%" mkdir "%resultsDir%"

:: Function to create LAMMPS input script
call :create_script >nul 2>&1

:: Run simulations in parallel
for %%s in (%seeds%) do (
    set "seedDir=%resultsDir%\seed_%%s"
    start "LAMMPS Seed %%s" cmd /c "cd /d "!seedDir!" && %lammpsExec% -in run_seed%%s.in > lammps_output_%%s.log 2>&1"
)

echo All simulations started in parallel...
echo Waiting for completion...

:: Wait for all processes to complete
:wait_loop
timeout /t 1 >nul
tasklist /fi "windowtitle eq LAMMPS Seed*" | find "cmd.exe" >nul
if not errorlevel 1 goto wait_loop

echo All simulations completed. Results are in %resultsDir%
pause
exit /b

:create_script
for %%s in (%seeds%) do (
    set "seedDir=%resultsDir%\seed_%%s"
    mkdir "!seedDir!" >nul 2>&1
    
    (
    echo # Initialization
    echo units           lj
    echo atom_style      full
    echo boundary        p p p
    echo read_data       %dataFile%
    echo.
    echo # Force field definitions 
    echo bond_style      fene
    echo bond_coeff      1 30.0 1.5 1.0 1.0 # k=30, r0=1.5, e=1, sigma=1
    echo.
    echo # Non-bonded interactions: LJ + Coulombic
    echo pair_style      lj/cut/coul/long 2.5 10.0
    echo pair_modify     shift yes
    echo.
    echo # Initial interaction parameters
    echo pair_coeff      * * 1.0 1.0
    echo.
    echo # Coulombic solver
    echo kspace_style    pppm 1.0e-4
    echo dielectric      1.0
    echo.
    echo special_bonds   lj 0.0 1.0 1.0 coul 0.0 1.0 1.0
    echo.
    echo # Neighbor lists
    echo neighbor        0.5 bin 
    echo neigh_modify    one 10000 every 1 delay 1 check yes
    echo.
    echo # System groups
    echo group           type_A type 1
    echo group           type_B type 2
    echo.
    echo # Energy minimization
    echo minimize        1.0e-4 1.0e-6 1000 10000
    echo.
    echo # Set up time integration
    echo timestep        0.005
    echo.
    echo # Temperature control with seed %%s
    echo variable        temp equal 2.0
    echo fix             1 all langevin 2.0 2.0 2.0 %%s
    echo fix             2 all nve
    echo.
    echo # Output settings
    echo thermo_style    custom step pe ke temp press
    echo thermo          1000
    echo.
    echo # ---- EQUILIBRATION RUN ----
    echo dump            eq_dump all custom 1000 "!seedDir!\equilibration_seed%%s.lammpstrj" id type q xu yu zu 
    echo dump_modify     eq_dump sort id 
    echo.
    echo run             50000
    echo undump          eq_dump
    echo reset_timestep 0
    echo.
    echo # Ready for production run 
    echo.
    echo #Compute RDFs between different types
    echo compute         rdf_all all rdf 100 1 2
    echo.
    echo # Compute gyration
    echo compute         rg_all all gyration
    echo.
    echo # MSD for dynamics
    echo compute         msd_all all msd com yes
    echo.
    echo # ---- OUTPUT FILES ----
    echo # RG output
    echo fix             rg_avg_all all ave/time 100 1 100 c_rg_all file "!seedDir!\rg_all_12_1_seed%%s.dat"
    echo # MSD for diffusion
    echo fix             msd_avg_all all ave/time 50 1 50 c_msd_all[4] file "!seedDir!\msd_all_12_1_seed%%s.dat" mode scalar
    echo.
    echo # ---- PRODUCTION RUN ----
    echo dump            prod_dump all custom 1000 "!seedDir!\production_seed%%s.lammpstrj" id type q x y z 
    echo dump_modify     prod_dump sort id
    echo.
    echo run             100000
    ) > "!seedDir!\run_seed%%s.in"
)
exit /b