import numpy as np
from baseclasses import AeroProblem
from reg_default_options import defaultFuncList

ap_tutorial_wing = AeroProblem(
    name="mdo_tutorial",
    alpha=1.8,
    mach=0.80,
    P=20000.0,
    T=220.0,
    areaRef=45.5,
    chordRef=3.25,
    beta=0.0,
    R=287.87,
    xRef=0.0,
    yRef=0.0,
    zRef=0.0,
    evalFuncs=defaultFuncList,
)

ap_tutorial_wing_laminar = AeroProblem(
    name="mdo_tutorial",
    alpha=1.8,
    beta=0.0,
    mach=0.50,
    P=137.0,
    T=293.15,
    R=287.87,
    areaRef=45.5,
    chordRef=3.25,
    xRef=0.0,
    yRef=0.0,
    zRef=0.0,
    evalFuncs=defaultFuncList,
)

ap_tutorial_wing_rotating = AeroProblem(
    name="mdo_tutorial",
    alpha=0.0,
    mach=0.1,
    rho=1.225,
    T=284.15,
    areaRef=45.5,
    chordRef=3.25,
    beta=0.0,
    R=287.87,
    xRef=0.0,
    yRef=0.0,
    zRef=0.0,
    evalFuncs=defaultFuncList,
)


ap_CRM = AeroProblem(
    name="CRM",
    alpha=1.8,
    mach=0.80,
    P=20000.0,
    T=220.0,
    areaRef=45.5,
    chordRef=3.25,
    beta=0.0,
    R=287.87,
    xRef=0.0,
    yRef=0.0,
    zRef=0.0,
    evalFuncs=defaultFuncList,
)


k = 0.0808
M = 0.6
gamma = 1.4
R = 287.085
T = 280.0
c = 1.0
alpha_m = 2.77  # 2.89 #2.77 #Modified numbers
alpha_0 = 2.34  # 2.41 #2.34

omega = 2 * M * np.sqrt(gamma * R * T) * k / c
deltaAlpha = -alpha_0 * np.pi / 180.0


ap_naca0012_time_acc = AeroProblem(
    name="0012pitching",
    alpha=alpha_m,
    mach=M,
    machRef=M,
    reynolds=4800000.0,
    reynoldsLength=c,
    T=T,
    R=R,
    areaRef=1.0,
    chordRef=c,
    evalFuncs=["cl", "cd", "cmz"],
    xRef=0.25,
    xRot=0.25,
    degreePol=0,
    coefPol=[0.0],
    degreeFourier=1,
    omegaFourier=omega,
    cosCoefFourier=[0.0, 0.0],
    sinCoefFourier=[deltaAlpha],
)

ap_naca64A010_time_spectral = AeroProblem(
    name="64A010pitchingTS",
    alpha=0.0,
    mach=0.796,
    reynolds=12.56e6,
    reynoldsLength=1.0,
    T=305.0,
    R=287.085,
    areaRef=1.0,
    chordRef=1.0,
    xRef=0.248,
    xRot=0.248,
    omegaFourier=112.59,
    evalFuncs=["cl", "cd", "cmz"],
)

ap_2D_conv_nozzle = AeroProblem(
    name="2D_conv_nozzle",
    alpha=00.0,
    mach=0.25,
    T=500,
    P=79326.7,
    areaRef=1.0,
    chordRef=2.0,
    R=287.87,
    evalFuncs=[
        "mdot",
        "mdot_up",
        "mdot_down",
        "mavgptot_up",
        "mavgptot_down",
        "aavgptot_up",
        "aavgptot_down",
        "mavgttot_up",
        "mavgttot_down",
        "mavgps_up",
        "mavgps_down",
        "aavgps_up",
        "aavgps_down",
        "mavgmn_up",
        "mavgmn_down",
        "thrust",
        "thrust_pressure",
        "thrust_viscous",
        "thrust_momentum",
    ],
)

ap_conic_conv_nozzle = AeroProblem(
    name="conic_conv_nozzle",
    alpha=90.0,
    mach=0.5,
    altitude=0.0,
    areaRef=1.0,
    chordRef=1.0,
    R=287.87,
    evalFuncs=[
        "mdot_up",
        "mdot_down",
        "mavgptot_up",
        "mavgptot_down",
        "aavgptot_up",
        "aavgptot_down",
        "mavgttot_up",
        "mavgttot_down",
        "mavgps_up",
        "mavgps_down",
        "aavgps_up",
        "aavgps_down",
    ],
)

ap_actuator_pipe = AeroProblem(
    name="actuator_pipe",
    alpha=00,
    mach=0.6,
    altitude=0.0,
    areaRef=1.0,
    chordRef=1.0,
    evalFuncs=[
        "mdot_in",
        "mdot_out",
        "aavgptot_in",
        "aavgptot_out",
        "mavgttot_in",
        "mavgttot_out",
        "aavgps_in",
        "aavgps_out",
        "area_in",
        "area_out",
        "mavgvx_in",
        "mavgvx_out",
        "forcexpressure_in",
        "forcexpressure_out",
        "forcexmomentum_in",
        "forcexmomentum_out",
        "flowpower_az",
    ],
)

ap_simple_cart_cube = AeroProblem(
    name="cube",
    V=32,  # m/s
    T=273 + 60,  # kelvin
    P=93e3,  # pa
    areaRef=1.0,  # m^2
    chordRef=1.0,  # m^2
    evalFuncs=["cd"],
    alpha=5.0,
    beta=5.0,
    xRef=0.0,
    yRef=0.0,
    zRef=0.0,
)
