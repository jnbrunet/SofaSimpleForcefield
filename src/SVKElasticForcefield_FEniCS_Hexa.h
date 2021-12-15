#pragma once

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseTopology.h>
#include <Eigen/Eigen>
#include "Hexahedron.h"


class SVKElasticForcefield_FEniCS_Hexa : public sofa::core::behavior::ForceField<sofa::defaulttype::Vec3Types> {
public:
    SOFA_CLASS(SVKElasticForcefield_FEniCS_Hexa, SOFA_TEMPLATE(sofa::core::behavior::ForceField, sofa::defaulttype::Vec3Types));

    // Aliases
    using Element = Hexahedron;
    using Coord = sofa::type::Vec3;

    template <class T>
    using Data = sofa::core::objectmodel::Data<T>;

    template <typename ObjectType>
    using Link = sofa::core::objectmodel::SingleLink<SVKElasticForcefield_FEniCS_Hexa, ObjectType, sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

    // Data structures

    struct GaussNode {
        Real weight = 0;
        Real jacobian_determinant = 0;
        Eigen::Matrix<double, Element::NumberOfNodes, 3> dN_dx = Eigen::Matrix<double, Element::NumberOfNodes, 3>::Zero();
    };

    // public methods
    SVKElasticForcefield_FEniCS_Hexa();

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
    /**
     * Get the complete tangent stiffness matrix as a compressed sparse matrix.
     *
     * \note This method will not reassembled the stiffness matrix. It will return
     *       the latest assembly done (usually during the latest Newton iteration).
     *       Use the update_stiffness() method to manually trigger a reassembly of
     *       the tangent stiffness matrix.
     * */
    auto K() const -> Eigen::SparseMatrix<Real> {
        using StorageIndex = typename Eigen::SparseMatrix<Real>::StorageIndex;

        // K is symmetric, so we only stored "one side" of the matrix.
        // But to accelerate the computation, coefficients were not
        // stored only in the upper or lower triangular part, but instead
        // in whatever triangular part (upper or lower) the first node
        // index of the element was. This means that a coefficient (i,j)
        // might be in the lower triangular part, while (k,l) is in the
        // upper triangular part. But no coefficient will be both in the
        // lower AND the upper part.

        std::vector<Eigen::Triplet<Real>> triplets;
        triplets.reserve(p_K.size()*2);
        for (StorageIndex k = 0; k < p_K.outerSize(); ++k) {
            for (typename Eigen::SparseMatrix<Real>::InnerIterator it(p_K, k); it; ++it) {
                const auto i = static_cast<StorageIndex>(it.row());
                const auto j = static_cast<StorageIndex>(it.col());
                const auto v = it.value();
                if (i != j) {
                    triplets.emplace_back(i, j, v);
                    triplets.emplace_back(j, i, v);
                } else {
                    triplets.emplace_back(i, i, v);
                }
            }
        }

        Eigen::SparseMatrix<Real> K;
        K.resize(p_K.rows(), p_K.cols());
        K.setFromTriplets(triplets.begin(), triplets.end());

        return K;
    }

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
    Eigen::SparseMatrix<Real> p_K;
};
