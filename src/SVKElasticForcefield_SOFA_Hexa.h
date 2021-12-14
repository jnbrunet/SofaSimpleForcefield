#pragma once

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseTopology.h>
#include <Eigen/Eigen>
#include "Hexahedron.h"


class SVKElasticForcefield_SOFA_Hexa : public sofa::core::behavior::ForceField<sofa::defaulttype::Vec3Types> {
public:
    SOFA_CLASS(SVKElasticForcefield_SOFA_Hexa, SOFA_TEMPLATE(sofa::core::behavior::ForceField, sofa::defaulttype::Vec3Types));

    // Aliases
    using Element = Hexahedron;
    using Coord = sofa::type::Vec3;

    template <class T>
    using Data = sofa::core::objectmodel::Data<T>;

    template <typename ObjectType>
    using Link = sofa::core::objectmodel::SingleLink<SVKElasticForcefield_SOFA_Hexa, ObjectType, sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

    // Data structures

    struct GaussNode {
        Real weight = 0;
        Real jacobian_determinant = 0;
        Eigen::Matrix<double, Element::NumberOfNodes, 3> dN_dx = Eigen::Matrix<double, Element::NumberOfNodes, 3>::Zero();
    };

    // public methods
    SVKElasticForcefield_SOFA_Hexa();

    /**
     * Initialize the forcefield by pre-computing the derivative of the shape function at each
     * Gauss point w.r.t the initial (undeformed) position of the mesh.
     */
    void init() override;

    /** StvK elasticity potential, i.e. W(x) = 1/2 lambda Tr^2(e) + mu e:e */
    double getPotentialEnergy (
            const sofa::core::MechanicalParams * mparams,
            const Data<sofa::type::vector<sofa::type::Vec3>> & d_x) const override;

    /** Linear elasticity residual, i.e. R(x) = dW(x)/dE */
    void addForce(
            const sofa::core::MechanicalParams * mparams,
            Data<sofa::type::vector<sofa::type::Vec3>> & d_f,
            const Data<sofa::type::vector<sofa::type::Vec3>>& d_x,
            const Data<sofa::type::vector<sofa::type::Vec3>>& d_v) override;

    /** Jacobian of the elastic residual, i.e. K(x) = dR(x)/dE */
    void addKToMatrix(
            sofa::defaulttype::BaseMatrix * matrix,
            double kFact,
            unsigned int & offset) override;

    /** Performs df = K*dx where K is the tangent stiffness matrix (see addKToMatrix) */
    void addDForce(
            const sofa::core::MechanicalParams* mparams,
            Data<sofa::type::vector<sofa::type::Vec3>> & d_df,
            const Data<sofa::type::vector<sofa::type::Vec3>> & d_dx) override;

    /** Draw the computational mesh as the simulation advance. Only for debug purposes. */
    void draw(const sofa::core::visual::VisualParams* vparams) override;

    /** Used to automatically define a minimal view box for the GUI */
    void computeBBox(const sofa::core::ExecParams* params, bool onlyVisible) override;

private:
    Data<double> d_youngModulus;
    Data<double> d_poissonRatio;
    Link<sofa::core::topology::BaseMeshTopology>   d_topology_container;          ///< Container of node indices per element
    std::vector<std::array<GaussNode, Element::NumberOfGaussNode>> p_gauss_nodes; ///< Set of Gauss nodes per elements
};
