#include "ink_simulation.h"
#include "point_agglomerate.h"

#include <iostream>
#include <cmath>
#include <random>
#include <limits>

// Helper function: generate a random number uniformly between a and b.
double random_uniform(double a, double b) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(a, b);
    return dis(gen);
}

// -----------------------
// InkSimulation Methods
// -----------------------

// Constructor: Initializes simulation parameters.
InkSimulation::InkSimulation(const std::vector<std::vector<Eigen::Vector3d>>& agglomeratePoints,
                             double timeStep,
                             double totalTime,
                             const Eigen::Vector3d& domainSize,
                             double particleRadius,
                             double particleMass,
                             const std::vector<double>& agglomerateDiameter,
                             double maxVelocity,
                             double maxAngularVelocity,
                             double targetPorosity,
                             bool animate)
    : BoundaryPhysics(domainSize),
      agglomeratePoints(agglomeratePoints),
      timeStep(timeStep),
      totalTime(totalTime),
      finalDomainSize(domainSize),
      domainSize(domainSize),
      particleRadius(particleRadius),
      particleMass(particleMass),
      agglomerateDiameter(agglomerateDiameter),
      maxVelocity(maxVelocity),
      maxAngularVelocity(maxAngularVelocity),
      targetPorosity(targetPorosity),
      animate(animate)
{}

// Destructor: Clean up allocated agglomerates.
InkSimulation::~InkSimulation() {
    for (auto agg : agglomerates) {
        delete agg;
    }
}

// Setup the simulation environment: Place agglomerates without overlap until target porosity is achieved.
void InkSimulation::setup() {
    std::cout << "Setting up simulation ..." << std::endl;

    std::vector<Eigen::Vector3d> positions;  // Stores the chosen COM positions for each agglomerate.
    usedAgglomerateIndices.clear();
    double totalSolidVolume = 0.0;

    // Compute volume per primary particle.
    double V_p = (4.0 / 3.0) * M_PI * std::pow(particleRadius, 3);
    std::vector<double> agglomerateVolumes;
    for (const auto& aggPts : agglomeratePoints) {
        int numParticles = aggPts.size();
        agglomerateVolumes.push_back(numParticles * V_p);
    }

    // Start with the desired final domain size as our initial domain.
    std::vector<double> currentDomainSize = { finalDomainSize(0), finalDomainSize(1), finalDomainSize(2) };
    double finalDomainVol = currentDomainSize[0] * currentDomainSize[1] * currentDomainSize[2];

    std::random_device rd;
    std::mt19937 gen(rd());

    // Continue adding agglomerates until the target porosity is achieved.
    while (true) {
        // Randomly select an agglomerate index.
        std::uniform_int_distribution<> disAgg(0, static_cast<int>(agglomeratePoints.size()) - 1);
        int aggIndex = disAgg(gen);
        const auto& aggPts = agglomeratePoints[aggIndex];
        double aggVolume = agglomerateVolumes[aggIndex];
        double aggDiam = agglomerateDiameter[aggIndex];

        bool placed = false;
        int attempt = 0;
        const int maxAttempts = 1000;

        while (!placed && attempt < maxAttempts) {
            // Generate a random position within the current domain.
            Eigen::Vector3d pos;
            pos(0) = random_uniform(0, currentDomainSize[0]);
            pos(1) = random_uniform(0, currentDomainSize[1]);
            pos(2) = random_uniform(0, currentDomainSize[2]);

            bool overlap = false;
            for (size_t j = 0; j < positions.size(); ++j) {
                Eigen::Vector3d otherPos = positions[j];
                Eigen::Vector3d separation = pos - otherPos;
                double distance = separation.norm();
                double minDistance = (aggDiam / 2.0) + (agglomerateDiameter[usedAgglomerateIndices[j]] / 2.0);
                if (distance < minDistance) {
                    overlap = true;
                    break;
                }
            }
            if (!overlap) {
                positions.push_back(pos);
                usedAgglomerateIndices.push_back(aggIndex);
                totalSolidVolume += aggVolume;
                placed = true;
            }
            ++attempt;
        }
        if (!placed) {
            // Expand domain size by 5% if placement failed.
            for (double &size : currentDomainSize) {
                size *= 1.05;
            }
            continue;
        }

        double solidFraction = totalSolidVolume / finalDomainVol;
        double porosity = 1.0 - solidFraction;
        std::cout << "\rCurrent porosity: " << porosity << " vs target: " << targetPorosity << std::flush;
        if (porosity <= targetPorosity) {
            break;
        }
    }
    std::cout << std::endl;

    // Store the initial (expanded) domain size for scaling.
    initialDomainSize = Eigen::Vector3d(currentDomainSize[0], currentDomainSize[1], currentDomainSize[2]);
    domainSize = initialDomainSize;

    // Now create agglomerate objects.
    for (size_t i = 0; i < positions.size(); ++i) {
        int idx = usedAgglomerateIndices[i];
        const auto& aggPts = agglomeratePoints[idx];

        // Random linear velocity vector.
        double linearSpeed = random_uniform(-maxVelocity, maxVelocity);
        Eigen::Vector3d linearDirection;
        linearDirection << random_uniform(-1,1), random_uniform(-1,1), random_uniform(-1,1);
        if(linearDirection.norm() < 1e-12)
            linearDirection = Eigen::Vector3d(1,0,0);
        linearDirection.normalize();
        Eigen::Vector3d initialVel = linearSpeed * linearDirection;

        // Random angular velocity vector.
        double angularSpeed = random_uniform(-maxAngularVelocity, maxAngularVelocity);
        Eigen::Vector3d angularDirection;
        angularDirection << random_uniform(-1,1), random_uniform(-1,1), random_uniform(-1,1);
        if(angularDirection.norm() < 1e-12)
            angularDirection = Eigen::Vector3d(1,0,0);
        angularDirection.normalize();
        Eigen::Vector3d initialOmega = angularSpeed * angularDirection;

        // Create a new PointAgglomerate using the selected agglomerate's particle positions.
        PointAgglomerate* agg = new PointAgglomerate(
            aggPts,
            particleMass,
            particleRadius,
            agglomerateDiameter[idx],
            positions[i],  // initial COM position
            initialVel,
            initialOmega,
            animate
        );
        agglomerates.push_back(agg);
    }
    std::cout << "Starting simulation ..." << std::endl;
}

// Adjust agglomerate positions according to the scaling factor (for domain reduction).
void InkSimulation::make_move(const Eigen::Vector3d& scalingFactor) {
    // Compute the center of the current domain.
    Eigen::Vector3d initialDomainCenter = domainSize / 2.0;
    // Scale the domain.
    domainSize = domainSize.cwiseProduct(scalingFactor);
    Eigen::Vector3d newDomainCenter = domainSize / 2.0;

    // Adjust each agglomerate's center-of-mass position.
    for (auto* agg : agglomerates) {
        agg->pos = (agg->pos - initialDomainCenter).cwiseProduct(scalingFactor) + newDomainCenter;
        // Update the global particle positions accordingly.
        agg->update_particle_positions();
    }
}

// detect_collisions: Check for collisions between agglomerates and resolve them.
// This simplified version uses bounding spheres first, then refines by checking per-particle distances.
void InkSimulation::detect_collisions() {
    int numAgg = static_cast<int>(agglomerates.size());
    for (int i = 0; i < numAgg; ++i) {
        for (int j = i + 1; j < numAgg; ++j) {
            PointAgglomerate* agg1 = agglomerates[i];
            PointAgglomerate* agg2 = agglomerates[j];

            // Use center-of-mass positions and bounding sphere diameters for broad-phase.
            double centerDistance = (agg1->pos - agg2->pos).norm();
            double minDist = 0.5 * (agg1->bsd + agg2->bsd);
            if (centerDistance < minDist) {
                // Narrow-phase: Check minimum distance between individual particles.
                double minParticleDistance = std::numeric_limits<double>::max();
                Eigen::Vector3d closestSep;
                int idx1 = -1, idx2 = -1;
                for (size_t m = 0; m < agg1->particle_positions.size(); ++m) {
                    for (size_t n = 0; n < agg2->particle_positions.size(); ++n) {
                        double d = (agg1->particle_positions[m] - agg2->particle_positions[n]).norm();
                        if (d < minParticleDistance) {
                            minParticleDistance = d;
                            closestSep = agg1->particle_positions[m] - agg2->particle_positions[n];
                            idx1 = static_cast<int>(m);
                            idx2 = static_cast<int>(n);
                        }
                    }
                }
                if (minParticleDistance < 2 * particleRadius) {
                    // Collision detected. Compute collision response.
                    Eigen::Vector3d r1_contact = agg1->particle_positions[idx1] - agg1->pos;
                    Eigen::Vector3d r2_contact = agg2->particle_positions[idx2] - agg2->pos;
                    Eigen::Vector3d n = closestSep.normalized();

                    // Velocities at the contact points.
                    Eigen::Vector3d v1 = (agg1->p / agg1->M) + agg1->omega.cross(r1_contact);
                    Eigen::Vector3d v2 = (agg2->p / agg2->M) + agg2->omega.cross(r2_contact);
                    Eigen::Vector3d v_rel = v1 - v2;
                    double v_rel_n = v_rel.dot(n);

                    if (v_rel_n < 0) {
                        // Compute impulse scalar.
                        Eigen::Vector3d r1_cross_n = r1_contact.cross(n);
                        Eigen::Vector3d r2_cross_n = r2_contact.cross(n);
                        double restitution = 1.0;
                        double J_num = -(1 + restitution) * v_rel_n;
                        double J_den = (1.0 / agg1->M + 1.0 / agg2->M +
                            r1_cross_n.dot(agg1->I_tensor.fullPivLu().solve(r1_cross_n)) +
                            r2_cross_n.dot(agg2->I_tensor.fullPivLu().solve(r2_cross_n)));
                        double J = J_num / J_den;
                        Eigen::Vector3d impulse = J * n;

                        // Update linear and angular momentum.
                        agg1->p += impulse;
                        agg2->p -= impulse;
                        agg1->L += r1_contact.cross(impulse);
                        agg2->L -= r2_contact.cross(impulse);

                        // Correct positions to resolve overlap.
                        double overlap = 2 * particleRadius - minParticleDistance;
                        Eigen::Vector3d correction = (overlap / 2.0) * n;
                        agg1->pos += correction * (agg2->M / (agg1->M + agg2->M));
                        agg2->pos -= correction * (agg1->M / (agg1->M + agg2->M));

                        // Update particle positions after correction.
                        agg1->update_particle_positions();
                        agg2->update_particle_positions();
                    }
                }
            }
        }
    }
}

// run(): Execute the simulation loop.
void InkSimulation::run() {
    // Setup environment and place agglomerates.
    setup();

    int totalSteps = static_cast<int>(totalTime / timeStep) + 1;
    // Set domainSize to the initial value.
    domainSize = initialDomainSize;

    double currentTime = 0.0;
    int step = 1;

    // Main simulation loop.
    // log domain size
    std::cout << "Initial domain size: " << domainSize.transpose() << std::endl;
    while (currentTime < totalTime) {
        // Scale the domain and adjust positions.
        // make_move(scalingFactor);
        // Detect and resolve collisions.
        detect_collisions();

        // Update each agglomerate's state.
        for (auto* agg : agglomerates) {
            apply_gravity_and_drag(agg, timeStep);
            agg->update(timeStep, domainSize);
            apply_boundary_conditions(agg);
        }


        if (animate) {
            write_positions_to_file(step);
        }

        currentTime += timeStep;
        step++;
        // std::cout << "\rSimulation time: " << currentTime << " / " << totalTime << std::flush;

        // log the position of the first agglomerate
        progress = currentTime / totalTime * 100;
        std::cout << "\rSimulation progress: " << progress << "%" << std::flush;
    }
    std::cout << std::endl;

    // Collect final results.
    positions.clear();
    rotationAngles.clear();
    for (auto* agg : agglomerates) {
        positions.push_back(agg->pos);
        rotationAngles.push_back(agg->get_cumulative_rotation_angles());
    }
}

// Accessors for simulation results.
const std::vector<Eigen::Vector3d>& InkSimulation::get_positions() const {
    return positions;
}

const std::vector<Eigen::Vector3d>& InkSimulation::get_rotation_angles() const {
    return rotationAngles;
}

const std::vector<int>& InkSimulation::get_used_agglomerate_indices() const {
    return usedAgglomerateIndices;
}
