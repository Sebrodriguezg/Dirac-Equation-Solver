import numpy as np
from .geometry import Grid
from .initial_state import GaussianPacket

class SimulationProblem:
    """
    Data class to hold all the parameters for a Dirac simulation.
    This object is constructed by the DiracProblemBuilder.
    """
    def __init__(self, grid, initial_state, potential, boundary_condition, time_step, total_time):
        self.grid = grid
        self.initial_state = initial_state
        self.potential = potential
        self.boundary_condition = boundary_condition
        self.time_step = time_step
        self.total_time = total_time

class DiracProblemBuilder:
    """
    Implements the Builder pattern to construct a SimulationProblem object.
    This provides a fluent and readable API for setting up a simulation.
    """
    def __init__(self):
        self._grid = None
        self._initial_state = None
        self._potential = None
        self._boundary_condition = None
        self._time_step = None
        self._total_time = None

    def set_grid(self, grid: Grid):
        """Sets the spatial grid for the simulation."""
        self._grid = grid
        return self

    def set_initial_state(self, initial_state: GaussianPacket):
        """Sets the initial spinor field configuration."""
        self._initial_state = initial_state
        return self

    def set_potential(self, potential):
        """Sets the scalar potential for the simulation."""
        self._potential = potential
        return self

    def set_boundary_condition(self, boundary_condition):
        """Sets the boundary condition strategy for the simulation."""
        self._boundary_condition = boundary_condition
        return self

    def set_time_parameters(self, time_step: float, total_time: float):
        """Sets the time step (dt) and total simulation time."""
        self._time_step = time_step
        self._total_time = total_time
        return self

    def build(self) -> SimulationProblem:
        """
        Constructs and returns the configured SimulationProblem object.
        Performs validation to ensure all required parameters are set.
        """
        if self._grid is None:
            raise ValueError("Grid must be set before building the problem.")
        if self._initial_state is None:
            raise ValueError("Initial state must be set before building.")
        if self._time_step is None or self._total_time is None:
            raise ValueError("Time parameters must be set.")
        if self._boundary_condition is None:
            raise ValueError("Boundary condition must be set.")

        # For now, potential is optional (defaults to free particle if not set)
        if self._potential is None:
            print("Warning: Potential not set. Defaulting to free particle (V=0).")

        return SimulationProblem(
            grid=self._grid,
            initial_state=self._initial_state,
            potential=self._potential,
            boundary_condition=self._boundary_condition,
            time_step=self._time_step,
            total_time=self._total_time,
        )

from . import _core, electron_mass

class DiracSolver:
    """
    The main solver class that orchestrates the simulation.
    It takes a SimulationProblem and initializes the C++ backend.
    """
    def __init__(self, problem: SimulationProblem):
        self.problem = problem

        # Generate the initial wave function on the grid
        psi_0 = problem.initial_state.evaluate_on_grid(problem.grid)

        # Create C++ Grid object
        cpp_grid = _core.Grid(problem.grid.shape, problem.grid.spacing)

        # Instantiate the C++ FDTD integrator
        self.integrator = _core.FDTDLeapfrogIntegrator(
            psi_0,
            cpp_grid,
            problem.potential,
            problem.boundary_condition,
            problem.time_step,
            electron_mass # Using the global constant
        )
        print(f"DiracSolver initialized with C++ '{self.integrator.get_name()}' engine.")


    def run_simulation(self):
        """
        Executes the time-evolution loop by calling the C++ step function.
        """
        num_steps = int(self.problem.total_time / self.problem.time_step)
        print(f"Running simulation for {num_steps} steps...")
        for i in range(num_steps):
            self.integrator.step()
            if (i + 1) % 10 == 0:
                print(f"  Step {i+1}/{num_steps} complete.")
        print("Simulation finished.")

    def get_psi(self):
        """Returns the current spinor field from the C++ backend."""
        return self.integrator.get_psi()

    def plot_probability_density(self):
        """
        Calculates and plots the probability density of the spinor field.
        """
        print("Plotting probability density...")