#! /usr/bin/python
# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)
# import funkcí z jiného adresáře
import os
import os.path

import pytest

path_to_script = os.path.dirname(os.path.abspath(__file__))
import unittest
import numpy as np
import sys

try:
    import skelet3d

    data3d = np.ones([3, 7, 9])
    data3d[:, 3, 3:6] = 0
    skelet3d.skelet3d(data3d)
    skelet3d_installed = True
    # skelet3d
except:
    skelet3d_installed = False
    logger.warning("skelet3d is not working")

try:
    import larcc

    larcc_installed = True
except:
    larcc_installed = False
    logger.warning("larcc is not working")

import teigen.tree
from teigen.tree import TreeBuilder

# There is some problem with VTK. Code seams to be fine but it fails
#  Generic Warning: In /tmp/vtk20150408-2435-1y7p97u/VTK-6.2.0/Common/Core/vtkObjectBase.cxx, line 93
#  Trying to delete object with non-zero reference count.
#  ERROR: In /tmp/vtk20150408-2435-1y7p97u/VTK-6.2.0/Common/Core/vtkObject.cxx, line 156
#  vtkObject (0x11a26e760): Trying to delete object with non-zero reference count.


VTK_MALLOC_PROBLEM = True


#

class TubeTreeTest(unittest.TestCase):
    def setUp(self):
        self.interactivetTest = False

    # interactivetTest = True

    @pytest.mark.LAR
    @unittest.skipIf(not ("larcc" in sys.modules), "larcc is not installed")
    def test_vessel_tree_lar(self):
        import teigen.tb_lar
        tvg = TreeBuilder(teigen.tb_lar.TBLar)
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        output = tvg.buildTree()  # noqa
        if self.interactiveTests:
            tvg.show()

    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_nothing(self):
        logger.debug("skelet3d_installed", skelet3d_installed)
        # import ipdb; ipdb.set_trace()
        self.assertTrue(False)

    # @unittest.skip("test debug")
    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_vessel_tree_vtk(self):
        tvg = TreeBuilder('vtk')
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        output = tvg.buildTree()  # noqa
        # tvg.show()
        # tvg.saveToFile("tree_output.vtk")

    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    @unittest.skipIf(not ("skelet3d" in sys.modules), "skelet3d is not installed")
    @unittest.skipIf(not skelet3d_installed, "skelet3d is not installed")
    def test_vessel_tree_vtk_from_skeleton(self):
        logger.debug("skelet3d_installed", skelet3d_installed)

        import skelet3d
        import skelet3d.skeleton_analyser
        import shutil

        fn_out = 'tree.vtk'
        if os.path.exists(fn_out):
            os.remove(fn_out)

        volume_data = np.zeros([3, 7, 9], dtype=np.int)
        volume_data[:, :, 1:3] = 1
        volume_data[:, 5, 2:9] = 1
        volume_data[:, 0:7, 5] = 1
        skelet = skelet3d.skelet3d(volume_data)

        skan = skelet3d.skeleton_analyser.SkeletonAnalyser(skelet, volume_data=volume_data, voxelsize_mm=[1, 1, 1])
        stats = skan.skeleton_analysis()

        tvg = TreeBuilder('vtk')
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        tvg.tree_data = stats
        output = tvg.buildTree()  # noqa
        tvg.saveToFile(fn_out)
        os.path.exists(fn_out)

    # TODO finish this test
    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_vessel_tree_vol(self):
        import teigen.tb_volume
        tvg = TreeBuilder(teigen.tb_volume.TBVolume)
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        output = tvg.buildTree()  # noqa
        # tvg.show()
        # if self.interactiveTests:
        #     tvg.show()

    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_import_new_vt_format(self):

        tvg = TreeBuilder()
        yaml_path = os.path.join(path_to_script, "vt_biodur.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [150, 150, 150]
        data3d = tvg.buildTree()

    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_cylinders_generator(self):
        from teigen.generators.cylinders import CylinderGenerator

        cg = CylinderGenerator()
        cg.run()

    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_vtk_tree(self):
        import numpy as np
        tree_data = {

        }
        element_number = 1
        area_size = 100
        radius = 5
        np.random.seed(0)
        pts1 = np.random.random([element_number, 3]) * (area_size - 4 * radius) + 2 * radius
        pts2 = np.random.random([element_number, 3]) * (area_size - 4 * radius) + 2 * radius
        for i in range(element_number):
            edge = {
                # "nodeA_ZYX_mm": vor3.vertices[simplex],
                # "nodeB_ZYX_mm": vor3.vertices[simplex],
                "nodeA_ZYX_mm": pts1[i],
                "nodeB_ZYX_mm": pts2[i],
                "radius_mm": radius
            }
            tree_data[i] = edge

        tvg = TreeBuilder('vtk')
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        # tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [area_size, area_size, area_size]
        tvg.tree_data = tree_data
        output = tvg.buildTree()  # noqa
        # tvg.show()
        tvg.saveToFile("test_tree_output.vtk")

    # @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_tree_generator(self):
        import numpy as np
        tree_data = {

        }
        element_number = 6
        np.random.seed(0)
        pts = np.random.random([element_number, 3]) * 100

        # construct voronoi
        import scipy.spatial
        import itertools
        vor3 = scipy.spatial.Voronoi(pts)

        # for i, two_points in enumerate(vor3.ridge_points):
        for i, simplex in enumerate(vor3.ridge_vertices):
            simplex = np.asarray(simplex)
            # fallowing line removes all ridges with oulayers
            simplex = simplex[simplex > 0]
            if np.all(simplex >= 0):

                x = vor3.vertices[simplex, 0]
                y = vor3.vertices[simplex, 1]
                z = vor3.vertices[simplex, 2]
                for two_points in itertools.combinations(simplex, 2):
                    edge = {
                        # "nodeA_ZYX_mm": vor3.vertices[simplex],
                        # "nodeB_ZYX_mm": vor3.vertices[simplex],
                        "nodeA_ZYX_mm": vor3.vertices[two_points[0]],
                        "nodeB_ZYX_mm": vor3.vertices[two_points[1]],
                        "radius_mm": 2
                    }
                    tree_data[i] = edge
            else:
                pass

        show_input_points = False
        if show_input_points:
            length = len(tree_data)
            for i in range(element_number):
                edge = {
                    #         #"nodeA_ZYX_mm": np.random.random(3) * 100,
                    "nodeA_ZYX_mm": pts[i - 1],
                    "nodeB_ZYX_mm": pts[i],
                    #         "nodeB_ZYX_mm": np.random.random(3) * 100,
                    "radius_mm": 1
                }
                tree_data[i + length] = edge

        tvg = TreeBuilder('vtk')
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        # tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        tvg.tree_data = tree_data
        output = tvg.buildTree()  # noqa
        # tvg.show()
        tvg.saveToFile("test_tree_output.vtk")

        tvgvol = TreeBuilder('vol')
        tvgvol.voxelsize_mm = [1, 1, 1]
        tvgvol.shape = [100, 100, 100]
        tvgvol.tree_data = tree_data
        outputvol = tvgvol.buildTree()
        tvgvol.saveToFile("tree_volume.pklz")
        # self.assertTrue(False)

    # @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_io3d(self):
        import io3d
        data3d = np.zeros([10, 10, 10])
        segmentation = np.zeros([10, 10, 10])

        data3d[2:7, :3:5, :6] = 100
        datap = {
            "data3d": data3d,
            # "segmentation": segmentation,
            "voxelsize_mm": [1, 1, 1]
        }
        io3d.write(datap, "file1.pklz")

    # @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_skimage_io_imsave(self):
        import skimage.io
        data3d = np.zeros([10, 10, 10])
        segmentation = np.zeros([10, 10, 10])

        data3d[2:7, :3:5, :6] = 100
        datap = {
            "data3d": data3d,
            # "segmentation": segmentation,
            "voxelsize_mm": [1, 1, 1]
        }
        skimage.io.imsave("skiamge.png", data3d[0])


if __name__ == "__main__":
    unittest.main()
