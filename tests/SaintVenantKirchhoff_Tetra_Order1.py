import sys
import unittest
from pathlib import Path

import Sofa
import meshio
import numpy as np
from scipy.sparse import csr_matrix

current_dir = Path(__file__).parent
site_packages_dir = (current_dir / '..' / '..' / 'lib' / 'python3' / 'site-packages').resolve()
sys.path.insert(0, str(site_packages_dir))
print(f'Adding {site_packages_dir} to sys.path')

beam_p1 = meshio.read(('../scenes/beam_p1.vtu'))
# Material
young_modulus = 3000
poisson_ratio = 0.3


def createScene(node):
    node.addObject('DefaultVisualManagerLoop')
    node.addObject('DefaultAnimationLoop')

    node.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
    node.addObject('RequiredPlugin',
                   pluginName="SofaOpenglVisual SofaSimpleForcefield SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine")

    node.addObject('RegularGridTopology', name="grid", min="-7.5 -7.5 0", max="7.5 7.5 80", n="9 9 21")
    node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-15",
                   relative_residual_tolerance_threshold="1e-10", printLog="1")
    node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")

    node.addObject('MechanicalObject', name="mo", src="@grid")
    node.addObject('TetrahedronSetTopologyContainer', name="topology")
    node.addObject('TetrahedronSetTopologyModifier')
    node.addObject('Hexa2TetraTopologicalMapping', input="@grid", output="@topology", swapping="1")

    node.addObject('SVKElasticForcefield_FEniCS', name="ff", youngModulus="3000", poissonRatio="0.3",
                   topology="@topology")

    node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
    node.addObject('FixedConstraint', indices="@fixed_roi.indices")
    node.addObject('BoxROI', name="top_roi", box="-7.5 -7.5 79.9 7.5 7.5 80.1")
    node.addObject('TriangleSetGeometryAlgorithms')
    node.addObject('TrianglePressureForceField', pressure="0 -10 0", topology="@topology",
                   triangleList="@top_roi.triangleIndices", showForces="1")


class TestHyperelasticForcefield():
    def assertMatrixEqual(self, A, B):
        return np.array_equal(A.todense(), B.todense())

    def assertMatrixNotEqual(self, A, B):
        return np.array_equal(A.todense(), B.todense())

    def test_stiffness_assembly(self):
        root = Sofa.Core.Node()
        createScene(root)
        Sofa.Simulation.init(root)
        x = np.array(root.mo.position.array(), dtype=np.float64, order='C', copy=False)
        K1 = csr_matrix(root.ff.K(), copy=True)

        # # Manually trigger the matrix assembly using a different position vector
        # x2 = x * 10
        # root.ff.assemble_stiffness(x2)
        # K2 = csr_matrix(root.ff.K(), copy=False)
        # with pytest.raises(TypeError):
        #     assert(self.assertMatrixNotEqual(K1, K2))

        # Manually trigger the matrix assembly using the same position vector
        root.ff.assemble_stiffness(x)
        K3 = csr_matrix(root.ff.K(), copy=False)
        assert (self.assertMatrixEqual(K1, K3))


if __name__ == '__main__':
    unittest.main()
