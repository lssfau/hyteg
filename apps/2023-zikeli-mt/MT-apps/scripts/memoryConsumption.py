import numpy as np

# configuration for Notebook and LSS
configs = {
        "2D-P1" : {"coarse": 1, "fine": 12, "err": 13},
        "3D-P1" : {"coarse": 1, "fine":  8, "err":  9},
        "2D-P2" : {"coarse": 1, "fine": 11, "err": 12},
        "3D-P2" : {"coarse": 1, "fine":  7, "err":  8},
        }
# configuration for Fritz
#configs = {
#        "2D-P1" : {"coarse": 1, "fine": 13, "err": 14},
#        "3D-P1" : {"coarse": 1, "fine":  9, "err": 10},
#        "2D-P2" : {"coarse": 1, "fine": 11, "err": 12},
#        "3D-P2" : {"coarse": 1, "fine":  7, "err":  8},
#        }
configs["2D-P1"]["DoFs"] =  np.array([
9.000000e+00,
2.500000e+01,
8.100000e+01,
2.890000e+02,
1.089000e+03,
4.225000e+03,
1.664100e+04,
6.604900e+04,
2.631690e+05,
1.050625e+06,
4.198401e+06,
1.678541e+07,
6.712525e+07,
2.684682e+08,
1.073807e+09,
4.295098e+09,
])
configs["3D-P1"]["DoFs"] =  np.array([
1.400000e+01,
6.300000e+01,
3.650000e+02,
2.457000e+03,
1.796900e+04,
1.373130e+05,
1.073345e+06,
8.487297e+06,
6.750285e+07,
5.384453e+08,
4.301262e+09,
])
configs["2D-P2"]["DoFs"] =  np.array([
2.500000e+01,
8.100000e+01,
2.890000e+02,
1.089000e+03,
4.225000e+03,
1.664100e+04,
6.604900e+04,
2.631690e+05,
1.050625e+06,
4.198401e+06,
1.678541e+07,
6.712525e+07,
2.684682e+08,
])
configs["3D-P2"]["DoFs"] =  np.array([
6.300000e+01,
3.650000e+02,
2.457000e+03,
1.796900e+04,
1.373130e+05,
1.073345e+06,
8.487297e+06,
6.750285e+07,
5.384453e+08,
])

functions = {
"general": { "minLvl": "coarse", "maxLvl": "fine"  , "precision": "bar" , "num": 3 },
"error"  : { "minLvl": "coarse", "maxLvl": "err"   , "precision": "real", "num": 2 },
"op-R"   : { "minLvl": "fine"  , "maxLvl": "fine"  , "precision": "bar" , "num": 1 },
"op-I"   : { "minLvl": "coarse", "maxLvl": "fine"  , "precision": "dot" , "num": 1 },
"init"   : { "minLvl": "coarse", "maxLvl": "fine"  , "precision": "real", "num": 2 },
"err-tmp": { "minLvl": "fine"  , "maxLvl": "err"   , "precision": "real", "num": 1 },
"IR-R"   : { "minLvl": "coarse", "maxLvl": "fine"  , "precision": "bar" , "num": 2 },
"IR-I"   : { "minLvl": "coarse", "maxLvl": "fine"  , "precision": "dot" , "num": 2 },
"GMG"    : { "minLvl": "coarse", "maxLvl": "fine"  , "precision": "dot" , "num": 1 },
"CG"     : { "minLvl": "coarse", "maxLvl": "coarse", "precision": "dot" , "num": 4 },
"Cheby"  : { "minLvl": "coarse", "maxLvl": "fine"  , "precision": "dot" , "num": 2 },
}

precisions = {
        "R64-S64": {"bar": 64, "dot": 64, "real": 64},
        "R64-S32": {"bar": 64, "dot": 32, "real": 64},
        "R64-S16": {"bar": 64, "dot": 16, "real": 64},
        "R32-S64": {"bar": 32, "dot": 64, "real": 64},
        "R32-S32": {"bar": 32, "dot": 32, "real": 64},
        "R32-S16": {"bar": 32, "dot": 16, "real": 64},
        "R16-S64": {"bar": 16, "dot": 64, "real": 64},
        "R16-S32": {"bar": 16, "dot": 32, "real": 64},
        "R16-S16": {"bar": 16, "dot": 16, "real": 64},
        }

solvers = {
        "IR": [
            "general",
            "error"  ,
            "op-R"   ,
            "op-I"   ,
#            "init"   ,
#            "err-tmp",
            "IR-R"   ,
            "IR-I"   ,
            "GMG"    ,
            "CG"     ,
            "Cheby"  ,
            ],
        "GMG": [
            "general",
            "error"  ,
            "op-R"   ,
            "op-I"   ,
#            "init"   ,
#            "err-tmp",
            "IR-R"   ,
            "IR-I"   ,
            "GMG"    ,
            "CG"     ,
            "Cheby"  ,
            ],
        }

interests = {
        " IR_R64-S32_3D-P1": { "solver": "IR" , "precision": "R64-S32", "config": "3D-P1", "mem": 0.0 },
        " IR_R32-S16_3D-P1": { "solver": "IR" , "precision": "R32-S16", "config": "3D-P1", "mem": 0.0 },
        "GMG_R64-S64_3D-P1": { "solver": "GMG", "precision": "R64-S64", "config": "3D-P1", "mem": 0.0 },
        "GMG_R32-S32_3D-P1": { "solver": "GMG", "precision": "R32-S32", "config": "3D-P1", "mem": 0.0 },
        "GMG_R16-S16_3D-P1": { "solver": "GMG", "precision": "R16-S16", "config": "3D-P1", "mem": 0.0 },
        " IR_R64-S32_2D-P1": { "solver": "IR" , "precision": "R64-S32", "config": "2D-P1", "mem": 0.0 },
        " IR_R32-S16_2D-P1": { "solver": "IR" , "precision": "R32-S16", "config": "2D-P1", "mem": 0.0 },
        "GMG_R64-S64_2D-P1": { "solver": "GMG", "precision": "R64-S64", "config": "2D-P1", "mem": 0.0 },
        "GMG_R32-S32_2D-P1": { "solver": "GMG", "precision": "R32-S32", "config": "2D-P1", "mem": 0.0 },
        "GMG_R16-S16_2D-P1": { "solver": "GMG", "precision": "R16-S16", "config": "2D-P1", "mem": 0.0 },
        }

for name, interest in interests.items():
    for fct in solvers[interest["solver"]]:
        minLvl = configs[interest["config"]][functions[fct]["minLvl"]]
        maxLvl = configs[interest["config"]][functions[fct]["maxLvl"]]
        for lvl in range(minLvl, maxLvl+1):
            interest["mem"] += configs[interest["config"]]["DoFs"][lvl] * precisions[interest["precision"]][functions[fct]["precision"]]

    #print(f"For {name} on max level {(configs[interest['config']]['fine']):2d} the memory consumption is:\t{(interest['mem']*1e-12):.5f} TB")
    print(f"For {name} on max level {(configs[interest['config']]['fine']):2d} the memory consumption is:\t{(interest['mem']*1e-9):.5f} GB")
