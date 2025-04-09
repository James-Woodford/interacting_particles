#ifndef POINT_AGGLOMERATE_H
#define POINT_AGGLOMERATE_H

#include <vector>
#include <Eigen/Dense>

class PointAgglomerate {
public:
    // Constructor: positions are the initial particle centers (in global coords)
    // The constructor computes the body frame (relative to the center of mass),
    // total mass, inertia tensor, etc.
    PointAgglomerate(const std::vector<Eigen::Vector3d>& positions,
                     double particle_mass,
                     double particle_radius,
                     double bounding_sphere_diameter,
                     const Eigen::Vector3d& initial_pos,
                     const Eigen::Vector3d& initial_vel,
                     const Eigen::Vector3d& initial_omega,
                     bool animate = false);

    // Destructor
    ~PointAgglomerate();

    // Update the agglomerate for a time step delta_t given the domain size (for wall collisions, if needed)
    void update(double delta_t, const Eigen::Vector3d& domain_size);

    // Update global particle positions: particle_positions = pos + r_body
    void update_particle_positions();

    // Recalculate the moment of inertia tensor from the current body frame positions.
    void recalculate_inertia_tensor();

    // Static: Convert an axis-angle pair to a rotation matrix (using Rodrigues’ formula)
    static Eigen::Matrix3d axis_angle_to_matrix(const Eigen::Vector3d& axis, double angle);

    // Return the cumulative rotation angles (Euler angles in XYZ order) extracted from the cumulative rotation matrix.
    Eigen::Vector3d get_cumulative_rotation_angles() const;

    // Compute total kinetic energy (translational + rotational)
    double get_kinetic_energy() const;

    // Return the current global particle positions.
    const std::vector<Eigen::Vector3d>& get_particle_positions() const;

    // Public members for simplicity (could be made private with getters/setters later)
    int num_particles;
    std::vector<Eigen::Vector3d> particle_positions; // Global positions for each particle
    std::vector<Eigen::Vector3d> r_body;             // Particle positions in body (local) frame relative to COM
    double m_p;            // Mass per particle
    double R;              // Particle radius
    double M;              // Total mass (num_particles * m_p)
    Eigen::Vector3d pos;   // Center-of-mass position
    Eigen::Vector3d p;     // Linear momentum
    Eigen::Vector3d L;     // Angular momentum
    Eigen::Vector3d omega; // Angular velocity
    Eigen::Matrix3d I_tensor; // Moment of inertia tensor
    Eigen::Matrix3d cumulative_rotation_matrix; // Accumulated rotation matrix

    double bsd;            // Bounding sphere diameter
    bool animate;          // Flag for animation (unused in core physics)

};

#endif // POINT_AGGLOMERATE_H
