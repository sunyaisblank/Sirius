@echo off
echo Starting Sirius System...
cd bin\Sirius.Build\bin\Release
if exist Sirius.exe (
    Sirius.exe
) else (
    echo Error: Sirius.exe not found!
    echo Please build the project first.
    pause
)
cd ..\..\..
