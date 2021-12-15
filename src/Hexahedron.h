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

        // Lagrangian interpolation for Linear Hexahedron elements
        // see: https://defelement.com/elements/lagrange.html 
        m <<    -u*v*w + u*v + u*w - u + v*w - v - w + 1, 
                u*(v*w - v - w + 1), 
                v*(u*w - u - w + 1), 
                u*v*(1 - w), 
                w*(u*v - u - v + 1), 
                u*w*(1-v), 
                v*w*(1-u),
                u*v*w;

        return m;
    }

    /// Derivatives of the shape function at each nodes w.r.t local coordinates and evaluated at (u, v, w)
    static inline auto dN(double u, double v, double w) -> Eigen::Matrix<double, NumberOfNodes, 3> {
        Eigen::Matrix<double, NumberOfNodes, 3> m;
        //      dL/du               dL/dv               dL/dw
        m <<    (v - 1)*(1 - w),   (w - 1)*(1 - u),    (u - 1)*(1 - v), 
                (v - 1)*(w - 1),    u*(w - 1),          u*(v - 1), 
                v*(w - 1),          (u - 1)*(w - 1),    v*(u - 1), 
                v*(1 - w),          u*(1 - w),          - u*v,
                w*(v - 1),          w*(u - 1),          (v - 1)*(u - 1),
                w*(1 - v),          -u*w,               u*(1 - v),
                - v*w,              w*(1 - u),          v*(1 - u),
                v*w,                u*w,                u*v;  
           
        return m;
    }

    /// Numerical integration points of the element
    static inline auto gauss_nodes() -> const std::array<GaussNode, NumberOfGaussNode> & {
        static const std::array<GaussNode, NumberOfGaussNode> gauss_nodes {
            GaussNode {Eigen::Vector3d(0.5 - 1/sqrt(3), 0.5 - 1/sqrt(3), 0.5 - 1/sqrt(3)), double(1)}, // Node 0
            GaussNode {Eigen::Vector3d(0.5 + 1/sqrt(3), 0.5 - 1/sqrt(3), 0.5 - 1/sqrt(3)), double(1)}, // Node 1
            GaussNode {Eigen::Vector3d(0.5 - 1/sqrt(3), 0.5 + 1/sqrt(3), 0.5 - 1/sqrt(3)), double(1)}, // Node 2
            GaussNode {Eigen::Vector3d(0.5 + 1/sqrt(3), 0.5 + 1/sqrt(3), 0.5 - 1/sqrt(3)), double(1)}, // Node 3
            GaussNode {Eigen::Vector3d(0.5 - 1/sqrt(3), 0.5 - 1/sqrt(3), 0.5 + 1/sqrt(3)), double(1)}, // Node 4
            GaussNode {Eigen::Vector3d(0.5 + 1/sqrt(3), 0.5 - 1/sqrt(3), 0.5 + 1/sqrt(3)), double(1)}, // Node 5
            GaussNode {Eigen::Vector3d(0.5 - 1/sqrt(3), 0.5 + 1/sqrt(3), 0.5 + 1/sqrt(3)), double(1)}, // Node 6
            GaussNode {Eigen::Vector3d(0.5 + 1/sqrt(3), 0.5 + 1/sqrt(3), 0.5 + 1/sqrt(3)), double(1)}  // Node 7

        };
        return gauss_nodes;
    }
};
