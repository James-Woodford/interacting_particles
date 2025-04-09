#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include "include/ink_simulation.h"
#include <iostream>
#include <vector>


namespace py = pybind11;

PYBIND11_MODULE(interacting_particles, m) {
    m.doc() = "Interacting Particles Simulation Module";

    py::class_<InkSimulation>(m, "InkSimulation")
        .def(py::init<
             const std::vector<std::vector<Eigen::Vector3d>>&,
             double,
             double,
             const Eigen::Vector3d&,
             double,
             double,
             const std::vector<double>&,
             double,
             double,
             double,
             bool
        >(),
        // kwargs:
             py::arg("agglomerate_points"),
             py::arg("time_step"),
             py::arg("total_time"),
             py::arg("domain_size"),
             py::arg("particle_radius"),
             py::arg("particle_mass"),
             py::arg("agglomerate_diameter"),
             py::arg("max_velocity"),
             py::arg("max_angular_velocity"),
             py::arg("target_porosity"),
             py::arg("animate"))

        // funcs
        .def("run", &InkSimulation::run)
        .def("get_positions", &InkSimulation::get_positions)
        .def("get_rotation_angles", &InkSimulation::get_rotation_angles)
        .def("get_used_agglomerate_indices", &InkSimulation::get_used_agglomerate_indices)
        ;
}
