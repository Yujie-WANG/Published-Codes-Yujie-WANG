# define the default global parameters
PARTS = Array[
    [602.0419, 1.48, 0.324324, 3E13],
    [2.5, 2.5, 1000.0, 0.0],
    [2.5, 2.5, 2000.0, 0.0],
    [2.5, 2.5, 2000.0 ,100.0, 1.67, 50.0, 1000.0, 50000.0, 0.05, 0.5, 10.0, 0.001, 1E6, 2E-3]
    ]
ENVIR = [25.0, 3.166-1.5, 40.0, 1000.0, 2.0, 740.0]
UP_P = 0.0
DROUGHT_P = 2.5
HISTORIES = [1, 1, 1]
LEGACIES = ones(3, 200, 2)
LEGACIES[:,:,1] = 0
ESFACTOR = 1.00

# arrays used for fitting issue
ARRAY_TAIR  = []
ARRAY_TLEAF = []
ARRAY_VPA   = []
ARRAY_CO2   = []
ARRAY_PAR   = []
ARRAY_WIND  = []
ARRAY_RABS  = []
ARRAY_PSOIL = []
ARRAY_AREA  = []
ARRAY_LABA  = []
ARRAY_VCMAX = []
ARRAY_JMAX  = []
ARRAY_PHOTO = []
ARRAY_TRANS = []
ARRAY_PLEAF = []

# values used in dynamic model
RAIN_MEAN = 0.01
RAIN_STDD = 0.001
RAIN_PROB = 0.10
VPD_MEAN = 1.5
VPD_STDD = 0.3
SOIL_VMAX = 0.1 # m3
SOIL_WMAX = 0.5 # ratio
TREE_AREA = 5E-4 #m2
DAY_LENGTH = 12.0
INI_PSOIL = 0.0
