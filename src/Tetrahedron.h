#pragma once

#include <Eigen/Eigen>
#include <array>

struct Tetrahedron {
    static constexpr auto NumberOfNodes = 4;
    static constexpr auto NumberOfGaussNode = 1;

    struct GaussNode {
        Eigen::Vector3d position;
        double weight;
    };

    /// Shape function of each nodes evaluated at local coordinates (u, v, w)
    static inline auto N(double u, double v, double w) -> Eigen::Matrix<double, NumberOfNodes, 1> {
        return {
                1 - u - v - w,  // Node 0
                u,              // Node 1
                v,              // Node 1
                w               // Node 3
        };
    }

    /// Derivatives of the shape function at each nodes w.r.t local coordinates and evaluated at (u, v, w)
    static inline auto dN(double /*u*/, double /*v*/, double /*w*/) -> Eigen::Matrix<double, NumberOfNodes, 3> {
        Eigen::Matrix<double, NumberOfNodes, 3> m;
        //  dL/du  dL/dv  dL/dw
        m <<   -1,    -1,    -1,   // Node 0
                1,     0,     0,   // Node 1
                0,     1,     0,   // Node 1
                0,     0,     1;   // Node 3
        return m;
    }

    /// Numerical integration points of the element
    static inline auto gauss_nodes() -> const std::array<GaussNode, NumberOfGaussNode> & {
        static const std::array<GaussNode, NumberOfGaussNode> gauss_nodes {
                GaussNode {Eigen::Vector3d(1/4., 1/4., 1/4.), double(1/6.)} // Node 0
        };
        return gauss_nodes;
    }
};