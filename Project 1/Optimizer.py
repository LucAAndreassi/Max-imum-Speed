from scipy.optimize import minimize
from Models import objective_function

# Initial guess for [speed, altitude] to start the search
initial_guess = [6, 6, 5, 5, 100, 100]  # Example initial guess for [wing_area, fuselage_length, SW2, SW3, V2, V3]

# Define ranges/bounds for variables: (min, max)
# Optimization parameters are: wing_area, fuselage_length, SW2, SW3, V2, V3
wing_area_bounds = (2, 10)
fuselage_length_bounds = (4, 7)
SW2_bounds = (1.0, 10)
SW3_bounds = (1.0, 10)
V2_bounds = (40, 120)
V3_bounds = (40, 120)

all_bounds = [
    wing_area_bounds,
    fuselage_length_bounds,
    SW2_bounds,
    SW3_bounds,
    V2_bounds,
    V3_bounds
]

# Run SciPy Nelder-Mead (excellent for simulation-based black-box problems)
optimizer_result = minimize(
    fun=objective_function, 
    x0=initial_guess, 
    bounds=all_bounds,
    method='Nelder-Mead'
)
if optimizer_result.success:
    best_wing_area, best_fuselage_length, best_SW2, best_SW3, best_V2, best_V3 = optimizer_result.x
    print("\n🎉 Optimization Successful!")
    print(f"Optimal Wing Area: {round(best_wing_area, 2)} ft²")
    print(f"Optimal Fuselage Length: {round(best_fuselage_length, 2)} ft")
    print(f"Optimal SW2: {round(best_SW2, 2)} lbs")
    print(f"Optimal SW3: {round(best_SW3, 2)} lbs")
    print(f"Optimal V2: {round(best_V2, 2)} m/s")
    print(f"Optimal V3: {round(best_V3, 2)} m/s")
    print(f"Lowest Achieved Cost Score: {round(optimizer_result.fun, 2)}")
else:
    print("\n❌ Optimizer failed to converge:", optimizer_result.message)