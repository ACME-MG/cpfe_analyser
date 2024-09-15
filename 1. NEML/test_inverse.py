
import math
from neml import elasticity, drivers
from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal
from neml.math import rotations

NUM_THREADS = 8
STRAIN_RATE = 1e-4
MAX_STRAIN  = 0.5

TAU_SAT  = 108.35
B        = 0.5840
TAU_0    = 120.21
GAMMA_0  = 3.3333e-5
N        = 2.5832

def deg_to_rad(degrees:float) -> float:
    if isinstance(degrees, list):
        return [deg_to_rad(d) for d in degrees]
    return degrees * math.pi / 180

def rad_to_deg(radians:float) -> float:
    if isinstance(radians, list):
        return [rad_to_deg(r) for r in radians]
    return radians * 180 / math.pi

def reorient_euler(euler:list) -> list:
    orientation = rotations.CrystalOrientation(euler[0], euler[1], euler[2], angle_type="radians", convention="bunge")
    inverse = orientation.inverse()
    new_euler = inverse.to_euler(angle_type="radians", convention="bunge")
    return list(new_euler)

lattice = crystallography.CubicLattice(1.0)
lattice.add_slip_system([1,1,0], [1,1,1])

eulers = [[0,135,45], [68,38,72], [110,75,80], [0,0,0]]
eulers = deg_to_rad(eulers)
orientations = [rotations.CrystalOrientation(*euler, angle_type="radians", convention="bunge") for euler in eulers]
weights = [1] * len(orientations)

e_model    = elasticity.CubicLinearElasticModel(205000, 138000, 126000, "components")
str_model  = slipharden.VoceSlipHardening(TAU_SAT, B, TAU_0)
slip_model = sliprules.PowerLawSlipRule(str_model, GAMMA_0, N)
i_model    = inelasticity.AsaroInelasticity(slip_model)
k_model    = kinematics.StandardKinematicModel(e_model, i_model)
sc_model   = singlecrystal.SingleCrystalModel(k_model, lattice, miter=16, max_divide=2, verbose=False)
pc_model   = polycrystal.TaylorModel(sc_model, orientations, nthreads=NUM_THREADS, weights=weights)
results    = drivers.uniaxial_test(pc_model, STRAIN_RATE, emax=MAX_STRAIN, nsteps=500, rtol=1e-6,
                                   atol=1e-10, miter=25, verbose=False, full_results=True)

orientations = pc_model.orientations(results["history"][0])
for orientation in orientations:
    euler = list(orientation.to_euler(angle_type="radians", convention="bunge"))
    inverse_euler = list(orientation.inverse().to_euler(angle_type="radians", convention="bunge"))
    print("=================")
    print("A:", rad_to_deg(euler))
    print("B:", rad_to_deg(inverse_euler))
    exit()
