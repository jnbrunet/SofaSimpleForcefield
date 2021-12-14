#pragma once

#include <Eigen/Eigen>
#include <array>

struct Hexahedron {
    static constexpr auto NumberOfNodes = 8;
    static constexpr auto NumberOfGaussNode = 8;

    struct GaussNode {
        Eigen::Vector3d position;
        double weight;
    };

    /// Shape function of each nodes evaluated at local coordinates (u, v, w)
    static inline auto N(double u, double v, double w) -> Eigen::Matrix<double, NumberOfNodes, 1> {
        Eigen::Matrix<double, NumberOfNodes, 1> m;
        m << (1/8.) * (1 - u) * (1 - v) * (1 - w),   // Node 0
             (1/8.) * (1 + u) * (1 - v) * (1 - w),   // Node 1
             (1/8.) * (1 + u) * (1 + v) * (1 - w),   // Node 2
             (1/8.) * (1 - u) * (1 + v) * (1 - w),   // Node 3
             (1/8.) * (1 - u) * (1 - v) * (1 + w),   // Node 4
             (1/8.) * (1 + u) * (1 - v) * (1 + w),   // Node 5
             (1/8.) * (1 + u) * (1 + v) * (1 + w),   // Node 6
             (1/8.) * (1 - u) * (1 + v) * (1 + w);   // Node 7
        return m;
    }

    /// Derivatives of the shape function at each nodes w.r.t local coordinates and evaluated at (u, v, w)
    static inline auto dN(double u, double v, double w) -> Eigen::Matrix<double, NumberOfNodes, 3> {
        Eigen::Matrix<double, NumberOfNodes, 3> m;
        //  dL/du  dL/dv  dL/dw
        m << -1/8. * (1 - v) * (1 - w),    -1/8. * (1 - u) * (1 - w),    -1/8. * (1 - u) * (1 - v),   // Node 0
             +1/8. * (1 - v) * (1 - w),    -1/8. * (1 + u) * (1 - w),    -1/8. * (1 + u) * (1 - v),   // Node 1
             +1/8. * (1 + v) * (1 - w),    +1/8. * (1 + u) * (1 - w),    -1/8. * (1 + u) * (1 + v),   // Node 2
             -1/8. * (1 + v) * (1 - w),    +1/8. * (1 - u) * (1 - w),    -1/8. * (1 - u) * (1 + v),   // Node 3
             -1/8. * (1 - v) * (1 + w),    -1/8. * (1 - u) * (1 + w),    +1/8. * (1 - u) * (1 - v),   // Node 4
             +1/8. * (1 - v) * (1 + w),    -1/8. * (1 + u) * (1 + w),    +1/8. * (1 + u) * (1 - v),   // Node 5
             +1/8. * (1 + v) * (1 + w),    +1/8. * (1 + u) * (1 + w),    +1/8. * (1 + u) * (1 + v),   // Node 6
             -1/8. * (1 + v) * (1 + w),    +1/8. * (1 - u) * (1 + w),    +1/8. * (1 - u) * (1 + v);   // Node 7
        return m;
    }

    /// Numerical integration points of the element
    static inline auto gauss_nodes() -> const std::array<GaussNode, NumberOfGaussNode> & {
        static const std::array<GaussNode, NumberOfGaussNode> gauss_nodes {
            GaussNode {Eigen::Vector3d(-1/sqrt(3), -1/sqrt(3), -1/sqrt(3)), double(1)}, // Node 0
            GaussNode {Eigen::Vector3d(+1/sqrt(3), -1/sqrt(3), -1/sqrt(3)), double(1)}, // Node 1
            GaussNode {Eigen::Vector3d(-1/sqrt(3), +1/sqrt(3), -1/sqrt(3)), double(1)}, // Node 2
            GaussNode {Eigen::Vector3d(+1/sqrt(3), +1/sqrt(3), -1/sqrt(3)), double(1)}, // Node 3
            GaussNode {Eigen::Vector3d(-1/sqrt(3), -1/sqrt(3), +1/sqrt(3)), double(1)}, // Node 4
            GaussNode {Eigen::Vector3d(+1/sqrt(3), -1/sqrt(3), +1/sqrt(3)), double(1)}, // Node 5
            GaussNode {Eigen::Vector3d(-1/sqrt(3), +1/sqrt(3), +1/sqrt(3)), double(1)}, // Node 6
            GaussNode {Eigen::Vector3d(+1/sqrt(3), +1/sqrt(3), +1/sqrt(3)), double(1)}  // Node 7
        };
        return gauss_nodes;
    }
};
