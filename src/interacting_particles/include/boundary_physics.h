// boundary_physics.h

#ifndef BOUNDARY_PHYSICS_H
#define BOUNDARY_PHYSICS_H

#include <array>
#include <Eigen/Dense>
#include <algorithm>
#include "point_agglomerate.h"

class BoundaryPhysics {
    public:
        BoundaryPhysics(const Eigen::Vector3d& domainSize,
                        const Eigen::Vector3d& gravity = Eigen::Vector3d(0, 0, -9000000.81),
                        double drag = 0.1,
                        double restitution_side = 0.3,
                        double restitution_bottom = 0.0)
            : domainSize(domainSize),
              gravity(gravity),
              dragCoefficient(drag),
              restitutionSide(restitution_side),
              restitutionBottom(restitution_bottom) {}

    protected:
        Eigen::Vector3d domainSize;
        Eigen::Vector3d gravity;
        double dragCoefficient;
        double restitutionSide;
        double restitutionBottom;

        void apply_gravity_and_drag(PointAgglomerate* agg, double dt) const {
            // Linear momentum: p = Mv => v = p / M
            Eigen::Vector3d velocity = agg->p / agg->M;

            // Apply gravity
            velocity[2] += gravity[2] * dt;


            // Apply drag (linear decay)
            velocity *= (1.0 - dragCoefficient * dt);

            // Update linear momentum
            agg->p = velocity * agg->M;
        }

        void apply_boundary_conditions(PointAgglomerate* agg) const {
            Eigen::Vector3d& pos = agg->pos;
            Eigen::Vector3d vel = agg->p / agg->M;

            // X boundaries
            if (pos[0] < 0 || pos[0] > domainSize[0]) {
                vel[0] *= -restitutionSide;
                pos[0] = std::clamp(pos[0], 0.0, domainSize[0]);
            }

            // Y boundaries
            if (pos[1] < 0 || pos[1] > domainSize[1]) {
                vel[1] *= -restitutionSide;
                pos[1] = std::clamp(pos[1], 0.0, domainSize[1]);
            }

            // Z bottom
            if (pos[2] < 0 ) {
                vel[2] *= -restitutionBottom;
                pos[2] = std::max(pos[2], 0.0);
            }
            // Z top, 1.0 restitution
            if (pos[2] > domainSize[2]) {
                vel[2] *= -1.0;
                pos[2] = std::min(pos[2], domainSize[2]);
            }

            // Update linear momentum
            agg->p = vel * agg->M;


    }
};

#endif // BOUNDARY_PHYSICS_H
