#pragma once

#include <Eigen/Eigen>
#include <array>

struct Tetrahedron {
    static constexpr auto NumberOfNodes = 10;
    static constexpr auto NumberOfGaussNode = 4;

    struct GaussNode {
        Eigen::Vector3d position;
        double weight;
    };

    /// Shape function of each nodes evaluated at local coordinates (u, v, w)
    static inline auto N(double u, double v, double w) -> Eigen::Matrix<double, NumberOfNodes, 1> {
        const auto l = 1 - u - v - w;
        Eigen::Matrix<double, NumberOfNodes, 1> m;
        m << l * (2*l - 1),   // Node 0
             u * (2*u - 1),   // Node 1
             v * (2*v - 1),   // Node 2
             w * (2*w - 1),   // Node 3
             4 * l * u,       // Node 4
             4 * u * v,       // Node 5
             4 * v * l,       // Node 6
             4 * l * w,       // Node 7
             4 * u * w,       // Node 8
             4 * v * w;       // Node 9

        return m;
    }

    /// Derivatives of the shape function at each nodes w.r.t local coordinates and evaluated at (u, v, w)
    static inline auto dN(double u, double v, double w) -> Eigen::Matrix<double, NumberOfNodes, 3> {
        const auto l = 1 - u - v - w;
        Eigen::Matrix<double, NumberOfNodes, 3> m;
        //       dL/du              dL/dv              dL/dw
        m <<   1 - (4 * l),       1 - (4 * l),       1 - (4 * l),   // Node 0
              (4 * u) - 1 ,          0       ,          0       ,   // Node 1
                    0     ,      (4 * v) - 1 ,          0       ,   // Node 2
                    0     ,          0       ,      (4 * w) - 1 ,   // Node 3
               4 * (l - u),      -4 * u      ,      -4 * u      ,   // Node 4
                4 * v     ,       4 * u      ,          0       ,   // Node 5
               -4 * v     ,       4 * (l - v),      -4 * v      ,   // Node 6
               -4 * w     ,      -4 * w      ,       4 * (l - w),   // Node 7
                4 * w     ,          0       ,       4 * u      ,   // Node 8
                    0     ,       4 * w      ,       4 * v      ;   // Node 9
        return m;
    }

    /// Numerical integration points of the element
    static inline auto gauss_nodes() -> const std::array<GaussNode, NumberOfGaussNode> & {
        static const std::array<GaussNode, NumberOfGaussNode> gauss_nodes {
                GaussNode {Eigen::Vector3d(0.1381966011250105, 0.1381966011250105, 0.1381966011250105), double(1/24.)}, // Node 0
                GaussNode {Eigen::Vector3d(0.1381966011250105, 0.1381966011250105, 0.5854101966249685), double(1/24.)}, // Node 1
                GaussNode {Eigen::Vector3d(0.1381966011250105, 0.5854101966249685, 0.1381966011250105), double(1/24.)}, // Node 2
                GaussNode {Eigen::Vector3d(0.5854101966249685, 0.1381966011250105, 0.1381966011250105), double(1/24.)}  // Node 3
        };
        return gauss_nodes;
    }
};
