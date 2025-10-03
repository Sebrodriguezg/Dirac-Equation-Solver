#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Required for vector bindings
#include "Grid.h"
#include "TimeIntegrator.h"
#include "potentials/Potential.h"
#include "potentials/FreeParticle.h"
#include "boundaries/BoundaryCondition.h"
#include "boundaries/PeriodicBoundary.h"

namespace py = pybind11;

// PYBIND11_MODULE defines the entry point for the Python extension module.
// The first argument is the name of the module (must match the one in CMakeLists.txt).
// The second argument, 'm', is a variable of type py::module_ which is the main interface
// for adding new bindings.
PYBIND11_MODULE(_core, m) {
    m.doc() = "pybind11 plugin for dirac_solver_core"; // optional module docstring

    // Bind the TimeIntegrator strategy interface
    // We create a trampoline class for pybind11 to handle virtual functions
    class PyTimeIntegrator : public TimeIntegrator {
    public:
        using TimeIntegrator::TimeIntegrator; // Inherit constructors
        void step() override { PYBIND11_OVERRIDE_PURE(void, TimeIntegrator, step); }
        std::string get_name() const override { PYBIND11_OVERRIDE_PURE(std::string, TimeIntegrator, get_name); }
    };

    py::class_<TimeIntegrator, PyTimeIntegrator, std::shared_ptr<TimeIntegrator>>(m, "TimeIntegrator")
        .def(py::init<>())
        .def("step", &TimeIntegrator::step)
        .def("get_name", &TimeIntegrator::get_name);

    // Bind the Grid class
    py::class_<Grid>(m, "Grid")
        .def(py::init<std::vector<size_t>, std::vector<double>>(),
             py::arg("shape"), py::arg("spacing"))
        .def("get_dim", &Grid::get_dim)
        .def("get_shape", &Grid::get_shape)
        .def("get_spacing", &Grid::get_spacing)
        .def("get_total_points", &Grid::get_total_points);

    // Bind the concrete FDTDLeapfrogIntegrator strategy
    py::class_<FDTDLeapfrogIntegrator, TimeIntegrator, std::shared_ptr<FDTDLeapfrogIntegrator>>(m, "FDTDLeapfrogIntegrator")
        .def(py::init<py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast>, const Grid&, std::shared_ptr<Potential>, std::shared_ptr<BoundaryCondition>, double, double>(),
             py::arg("initial_psi"), py::arg("grid"), py::arg("potential"), py::arg("boundary_condition"), py::arg("dt"), py::arg("mass"))
        .def("step", &FDTDLeapfrogIntegrator::step)
        .def("get_psi", &FDTDLeapfrogIntegrator::get_psi, "Returns the current spinor field as a NumPy array");

    // Bind the BoundaryCondition base class (Strategy)
    class PyBoundaryCondition : public BoundaryCondition {
    public:
        using BoundaryCondition::BoundaryCondition; // Inherit constructors
        void apply(std::vector<Dirac::Spinor>& psi) const override {
            PYBIND11_OVERRIDE_PURE(void, BoundaryCondition, apply, psi);
        }
        std::string get_name() const override {
            PYBIND11_OVERRIDE_PURE(std::string, BoundaryCondition, get_name);
        }
    };

    py::class_<BoundaryCondition, PyBoundaryCondition, std::shared_ptr<BoundaryCondition>>(m, "BoundaryCondition")
        .def(py::init<>())
        .def("apply", &BoundaryCondition::apply, "Applies the boundary condition")
        .def("get_name", &BoundaryCondition::get_name);

    // Bind the concrete PeriodicBoundary condition
    py::class_<PeriodicBoundary, BoundaryCondition, std::shared_ptr<PeriodicBoundary>>(m, "PeriodicBoundary")
        .def(py::init<>())
        .def("get_name", &PeriodicBoundary::get_name);

    // Bind the Potential base class (Template)
    class PyPotential : public Potential {
    public:
        using Potential::Potential; // Inherit constructors
        double evaluate(const std::vector<double>& r) const override {
            PYBIND11_OVERRIDE_PURE(double, Potential, evaluate, r);
        }
        std::string get_name() const override {
            PYBIND11_OVERRIDE_PURE(std::string, Potential, get_name);
        }
    };

    py::class_<Potential, PyPotential, std::shared_ptr<Potential>>(m, "Potential")
        .def(py::init<>())
        .def("evaluate", &Potential::evaluate, "Evaluates the potential at a given coordinate", py::arg("r"))
        .def("get_name", &Potential::get_name);

    // Bind the concrete FreeParticle potential
    py::class_<FreeParticle, Potential, std::shared_ptr<FreeParticle>>(m, "FreeParticle")
        .def(py::init<>())
        .def("get_name", &FreeParticle::get_name); // evaluate() is inherited and works
}