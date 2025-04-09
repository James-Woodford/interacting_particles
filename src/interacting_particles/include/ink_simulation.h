#ifndef INK_SIMULATION_H
#define INK_SIMULATION_H

#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <sstream>
#include <iostream>
#include "point_agglomerate.h"
#include "boundary_physics.h"

class InkSimulation : public BoundaryPhysics {
    // This class simulates the dynamics of agglomerates in a fluid medium.
    // It handles the initialization, simulation loop, and collision detection
    // for a set of agglomerates represented by their primary particles.
public:
    // Constructor: initialize simulation parameters.
    //   agglomeratePoints: Each agglomerate's initial particle positions (global coordinates).
    //   timeStep: Simulation time step.
    //   totalTime: Total simulation duration.
    //   domainSize: Final desired domain dimensions.
    //   particleRadius: Radius of each spherical primary particle.
    //   particleMass: Mass of each primary particle.
    //   agglomerateDiameter: Effective diameters of the agglomerates.
    //   maxVelocity: Maximum initial linear speed.
    //   maxAngularVelocity: Maximum initial angular speed.
    //   targetPorosity: Target void fraction.
    //   animate: Flag to indicate if animation-related code should run.
    InkSimulation(const std::vector<std::vector<Eigen::Vector3d>>& agglomeratePoints,
                  double timeStep,
                  double totalTime,
                  const Eigen::Vector3d& domainSize,
                  double particleRadius,
                  double particleMass,
                  const std::vector<double>& agglomerateDiameter,
                  double maxVelocity,
                  double maxAngularVelocity,
                  double targetPorosity,
                  bool animate = false);

    // Destructor: cleans up dynamically allocated agglomerates.
    ~InkSimulation();

    // Run the simulation.
    // This method sets up the initial state, iterates the simulation loop,
    // and computes final results (positions, rotation angles, used agglomerate indices).
    void run();

    // Accessor methods for simulation results.
    const std::vector<Eigen::Vector3d>& get_positions() const;
    const std::vector<Eigen::Vector3d>& get_rotation_angles() const;
    const std::vector<int>& get_used_agglomerate_indices() const;

    // Detect and resolve collisions between agglomerate particles.
    void detect_collisions();

    // Setup the simulation environment: place agglomerates without overlap until target porosity is reached.
    void setup();

    // Adjust agglomerate positions based on a scaling factor (for domain size reduction).
    // The scaling factor is provided as an Eigen::Vector3d (one factor per axis).
    void make_move(const Eigen::Vector3d& scalingFactor);

    std::vector<Eigen::Vector3d> getAgglomeratePositions()
    {
        std::vector<Eigen::Vector3d> currentPositions;
        for (auto* agg : agglomerates) {
            currentPositions.push_back(agg->pos);
        }
        return currentPositions;

    }

    void write_positions_to_file(int timestep)
    {
        // Create output folder if it doesn't exist
        std::filesystem::create_directories("positions_files");

        // Filename with padded timestep number
        std::ostringstream filename;
        filename << "positions_files/positions_"
                 << std::setfill('0') << std::setw(4) << timestep
                 << ".txt";

        std::ofstream file(filename.str());
        if (!file.is_open()) {
            std::cout << "Error: could not open file for writing: " << filename.str() << std::endl;
            return;
        }

        // Write positions: one line per particle: x y z
        for (const auto& pos : getAgglomeratePositions()) {
            file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
        }

        file.close();
    }

private:
    // Simulation parameters.
    std::vector<std::vector<Eigen::Vector3d>> agglomeratePoints;  // Each agglomerate's list of primary particle positions.
    std::vector<double> agglomerateDiameter;                      // Effective diameters of the agglomerates.
    double timeStep;
    double totalTime;
    Eigen::Vector3d finalDomainSize;   // Desired final domain dimensions.
    Eigen::Vector3d domainSize;        // Current domain dimensions.
    double particleRadius;
    double particleMass;
    double maxVelocity;
    double maxAngularVelocity;
    double targetPorosity;
    bool animate;

    // Simulation state.
    std::vector<PointAgglomerate*> agglomerates;  // Collection of created agglomerates.
    std::vector<int> usedAgglomerateIndices;      // Indices of agglomerates that were placed.
    Eigen::Vector3d initialDomainSize;            // Initial (possibly expanded) domain size for scaling.

    // Results (collected after simulation run).
    std::vector<Eigen::Vector3d> positions;         // Final center-of-mass positions for each agglomerate.
    std::vector<Eigen::Vector3d> rotationAngles;    // Cumulative rotation angles (e.g., Euler angles) for each agglomerate.






};

#endif // INK_SIMULATION_H