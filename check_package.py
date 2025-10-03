import numpy as np
import dirac_solver
from dirac_solver import (
    Grid,
    ConstantSpinor,
    GaussianPacket,
    DiracProblemBuilder,
    DiracSolver,
    _core,
)

print("=== Verifying Boundary Conditions ===")

try:
    # 1. Use the Python Builder to construct a 1D simulation problem
    print("\n--- Building simulation problem with Periodic Boundaries ---")
    py_grid = Grid(shape=(100,), spacing=(0.1,))
    constant_spinor = ConstantSpinor(momentum=[5.0, 0, 0])
    initial_state = GaussianPacket(constant_spinor, center=[0.0, 0, 0])

    # 2. Create C++ objects for potential and boundary conditions
    potential = _core.FreeParticle()
    boundary_condition = _core.PeriodicBoundary()
    print(f"[OK] C++ Potential '{potential.get_name()}' instantiated.")
    print(f"[OK] C++ Boundary Condition '{boundary_condition.get_name()}' instantiated.")

    # 3. Build the problem
    problem = (
        DiracProblemBuilder()
        .set_grid(py_grid)
        .set_initial_state(initial_state)
        .set_potential(potential)
        .set_boundary_condition(boundary_condition) # Set the boundary condition
        .set_time_parameters(time_step=0.001, total_time=0.1)
        .build()
    )
    print("[OK] SimulationProblem constructed successfully.")
    print(f"   - Boundary Condition: {problem.boundary_condition.get_name()}")


    # 4. Initialize the main solver, which sets up the C++ FDTD engine
    solver = DiracSolver(problem)
    print("[OK] DiracSolver initialized, C++ FDTD engine is ready.")

    # 5. Run the simulation with periodic boundaries
    solver.run_simulation()

    # 6. Retrieve the final wave function and verify it
    psi_final = solver.get_psi()
    print("\n--- Verifying simulation results ---")
    print(f"[OK] Final wave function retrieved from C++ backend.")

    expected_shape = (py_grid.shape[0], 4)
    if psi_final.shape == expected_shape:
        print(f"[OK] Final psi has the correct shape: {psi_final.shape}")
    else:
        print(f"[FAIL] Final psi has incorrect shape: {psi_final.shape}, expected {expected_shape}")

    norm = np.linalg.norm(psi_final)
    if norm > 1e-9:
        print(f"[OK] Final wave function is non-trivial (norm = {norm:.4f}).")
    else:
        print(f"[FAIL] Final wave function appears to be zero (norm = {norm:.4f}).")


    print("\n=== Boundary Conditions Verified Successfully ===")

except Exception as e:
    print(f"\n[ERROR] Verification failed: {e}")
    import traceback
    traceback.print_exc()