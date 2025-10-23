# Before trying to run MC_AdvectionDiffusionConvTest.cpp or MC_PureAdvectionConvTest.cpp:
# Run this file and run CMake configure afterwards to create links to the parameter files in the build folder.

timeLevels = [0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125, 0.0015625, 0.00078125, 0.000390625, 0.0001953125, 9.765625e-05, 4.8828125e-05, 2.44140625e-05]

def createParameterString(_gridLevel: int, _dt: float):
    return (
    f'Parameters\n'
    f'{{\n'
    f'   minLevel {_gridLevel};\n'
    f'   maxLevel {_gridLevel};\n'
    f'   vtk false;\n'
    f'   rel_tolerance 1e-10;\n'
    f'   abs_tolerance 1e-12;\n'
    f'   t_init 3.5;\n'
    f'   run_time 1.0;\n'
    f'   delta_t {_dt};\n'
    f'   k 3;   \n'
    f'   \n'
    f'   useExtrapolation true;\n'
    f'   VelocityExtrapolationOrder 1;\n'
    f'   TemperatureExtrapolationOrder 1;\n'
    f'   \n'
    f'   BDFOrder 2;\n'
    f'   CrankNicolson false;\n'
    f'   MMOC true;\n'
    f'}}'
    )
    
for gridLevel in range(8):
    for timeLevel in range(13):
        dt = timeLevels[timeLevel]
        with open(f"BDF2_Level{gridLevel}_TimeLevel{timeLevel}.prm", "w") as f_out:
            f_out.write(createParameterString(gridLevel, dt))
