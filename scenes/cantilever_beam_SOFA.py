# Required import for python
import Sofa, SofaCaribou
import meshio
import numpy as np


# Choose in your script to activate or not the GUI
USE_GUI = True


# Import beam P2 elements
beam_p2 = meshio.read(('./beam_p2.vtu'))
# Material
young_modulus = 3000
poisson_ratio = 0.3


def main():
    import SofaRuntime
    import Sofa.Gui
    SofaRuntime.importPlugin("SofaOpenglVisual")
    SofaRuntime.importPlugin("SofaImplicitOdeSolver")
    SofaRuntime.importPlugin("SofaLoader")

    root = Sofa.Core.Node("root")
    createScene(root)
    Sofa.Simulation.init(root)

    if not USE_GUI:
        for iteration in range(10):
            Sofa.Simulation.animate(root, root.dt.value)
    else:
        Sofa.Gui.GUIManager.Init("myscene", "qglviewer")
        Sofa.Gui.GUIManager.createGUI(root, __file__)
        Sofa.Gui.GUIManager.SetDimension(1080, 1080)
        Sofa.Gui.GUIManager.MainLoop(root)
        Sofa.Gui.GUIManager.closeGUI()


def createScene(root):
    root.addObject('DefaultVisualManagerLoop')
    root.addObject('DefaultAnimationLoop')

    root.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
    root.addObject('RequiredPlugin', pluginName="SofaOpenglVisual SofaSimpleForcefield SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine")
    root.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-15", relative_residual_tolerance_threshold="1e-10", printLog="1")
    root.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")

    root.addChild("tet10")
    root.tet10.addObject('MechanicalObject', name='mo', position=beam_p2.points.tolist(), showObject=True, showObjectScale=5)
    root.tet10.addObject('CaribouTopology', name='volumetric_topology', template='Tetrahedron10', indices=beam_p2.cells_dict['tetra10'].tolist())
    root.tet10.addObject('SaintVenantKirchhoffMaterial', young_modulus=young_modulus, poisson_ratio=poisson_ratio)
    root.tet10.addObject('HyperelasticForcefield', printLog=True)
    # root.tet10.addObject('SVKElasticForcefield_SOFA_Tetra_Order2', youngModulus=young_modulus, poissonRatio=poisson_ratio, topology="@volumetric_topology")
    root.tet10.addObject('BoxROI', name='fixed_roi', box=[-7.5, -7.5, -0.9, 7.5, 7.5, 0.1])
    root.tet10.addObject('FixedConstraint', indices='@fixed_roi.indices')
    root.tet10.addObject('CaribouTopology', name='surface_topology', template='Triangle6', indices=beam_p2.cells_dict['triangle6'][np.array(np.ma.masked_equal(beam_p2.cell_data['gmsh:physical'][0], 2).mask)].tolist())
    root.tet10.addObject('TractionForcefield', traction=[0, -10, 0], slope=1, topology='@surface_topology', printLog=True)

    return root


# Function used only if this script is called from a python environment
if __name__ == '__main__':
    main()