//
// Created by james on 05/04/25.
//

#include "point_agglomerate.h"
#include <cmath>
#include <iostream>

// Constructor
PointAgglomerate::PointAgglomerate(const std::vector<Eigen::Vector3d>& positions,
                                   double particle_mass,
                                   double particle_radius,
                                   double bounding_sphere_diameter,
                                   const Eigen::Vector3d& initial_pos,
                                   const Eigen::Vector3d& initial_vel,
                                   const Eigen::Vector3d& initial_omega,
                                   bool animate)
    : m_p(particle_mass),
      R(particle_radius),
      bsd(bounding_sphere_diameter),
      pos(initial_pos),
      animate(animate),
      cumulative_rotation_matrix(Eigen::Matrix3d::Identity())
{
    // Number of particles
    num_particles = positions.size();

    // Total mass
    M = num_particles * m_p;

    // Linear momentum
    p = M * initial_vel;

    // Angular momentum initialized to zero; will update after computing inertia tensor
    L = Eigen::Vector3d::Zero();

    // Set initial angular velocity
    omega = initial_omega;

    // Store original particle positions (in global coordinates)
    // We will compute the body frame positions relative to the COM computed from these.
    std::vector<Eigen::Vector3d> initial_particle_positions = positions;

    // Compute the center-of-mass (body frame) from given positions
    Eigen::Vector3d r_cm_body = Eigen::Vector3d::Zero();
    for (const auto& r : initial_particle_positions) {
        r_cm_body += r;
    }
    r_cm_body /= static_cast<double>(num_particles);

    // Compute body frame positions: each particle relative to the computed COM.
    r_body.resize(num_particles);
    for (size_t i = 0; i < num_particles; ++i) {
        r_body[i] = initial_particle_positions[i] - r_cm_body;
    }

    // Set the global particle positions based on the provided initial COM (pos)
    // This effectively "re-centers" the agglomerate at the given initial position.
    particle_positions.resize(num_particles);
    for (size_t i = 0; i < num_particles; ++i) {
        particle_positions[i] = pos + r_body[i];
    }

    // Compute the inertia tensor: I_tensor += m_p*(|r|^2*I - r*r^T)
    I_tensor = Eigen::Matrix3d::Zero();
    for (const auto& r : r_body) {
        I_tensor += m_p * (r.dot(r) * Eigen::Matrix3d::Identity() - r * r.transpose());
    }

    // Initialize angular momentum from initial angular velocity:
    L = I_tensor * omega;
}

// Destructor
PointAgglomerate::~PointAgglomerate() {
    // Nothing to free since we use STL containers and Eigen objects.
}

// Update the agglomerate for a time step delta_t (includes translation, rotation, and wall collision handling if desired)
void PointAgglomerate::update(double delta_t, const Eigen::Vector3d& domain_size) {
    // Update COM position based on current momentum (v_cm = p/M)
    Eigen::Vector3d v_cm = p / M;
    pos += v_cm * delta_t;

    // Update angular velocity based on current angular momentum.
    // If I_tensor is invertible, solve I * omega = L.
    if (std::abs(I_tensor.determinant()) > 1e-12) {
        omega = I_tensor.fullPivLu().solve(L);
    } else {
        omega = Eigen::Vector3d::Zero();
    }

    // Rotation
    double omega_norm = omega.norm();
    if (omega_norm > 1e-12) {
        double angle = omega_norm * delta_t;
        Eigen::Vector3d axis = omega / omega_norm;

        Eigen::Matrix3d delta_R = axis_angle_to_matrix(axis, angle);
        for (size_t i = 0; i < r_body.size(); ++i) {
            r_body[i] = delta_R * r_body[i];
        }
        cumulative_rotation_matrix = delta_R * cumulative_rotation_matrix;
    }

    update_particle_positions();

    // periodic boundary conditions tba
}

// Updates global particle positions: global = COM + relative (body-frame) position.
void PointAgglomerate::update_particle_positions() {
    for (size_t i = 0; i < num_particles; ++i) {
        particle_positions[i] = pos + r_body[i];
    }
}

// Recalculate the inertia tensor based on current body frame positions.
void PointAgglomerate::recalculate_inertia_tensor() {
    I_tensor = Eigen::Matrix3d::Zero();
    for (const auto& r : r_body) {
        I_tensor += m_p * (r.dot(r) * Eigen::Matrix3d::Identity() - r * r.transpose());
    }
}

// Static: Convert an axis-angle pair to a rotation matrix
Eigen::Matrix3d PointAgglomerate::axis_angle_to_matrix(const Eigen::Vector3d& axis, double angle) {
    Eigen::Vector3d n = axis.normalized();
    double cos_angle = std::cos(angle);
    double sin_angle = std::sin(angle);
    double one_minus_cos = 1.0 - cos_angle;
    double x = n(0), y = n(1), z = n(2);

    Eigen::Matrix3d R;
    R << cos_angle + x*x*one_minus_cos,      x*y*one_minus_cos - z*sin_angle,  x*z*one_minus_cos + y*sin_angle,
         y*x*one_minus_cos + z*sin_angle,      cos_angle + y*y*one_minus_cos,    y*z*one_minus_cos - x*sin_angle,
         z*x*one_minus_cos - y*sin_angle,      z*y*one_minus_cos + x*sin_angle,  cos_angle + z*z*one_minus_cos;
    return R;
}

// Return Euler angles (XYZ convention) from the cumulative rotation matrix.
Eigen::Vector3d PointAgglomerate::get_cumulative_rotation_angles() const {
    // Eigen's eulerAngles() returns angles in radians.
    // The following extracts intrinsic rotations about X, Y, Z in that order.
    return cumulative_rotation_matrix.eulerAngles(0, 1, 2);
}

// Compute the total kinetic energy (translational + rotational)
double PointAgglomerate::get_kinetic_energy() const {
    Eigen::Vector3d v_cm = p / M;
    double KE_trans = 0.5 * M * v_cm.squaredNorm();
    double KE_rot = 0.5 * omega.dot(I_tensor * omega);
    return KE_trans + KE_rot;
}

// Return the current global particle positions.
const std::vector<Eigen::Vector3d>& PointAgglomerate::get_particle_positions() const {
    return particle_positions;
}


