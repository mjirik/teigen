#! /usr/bin/python
# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)
# import funkcí z jiného adresáře
import os
import os.path

from nose.plugins.attrib import attr
path_to_script = os.path.dirname(os.path.abspath(__file__))
import unittest
import numpy as np
import sys

try:
    import skelet3d
    data3d = np.ones([3,7,9])
    data3d[:,3,3:6] = 0
    skelet3d.skelet3d(data3d)
    # skelet3d
except:
    pass
try:
    import larcc
except:
    pass
import teigen.tree
from teigen.tree import TreeBuilder

#

class TubeTreeTest(unittest.TestCase):
    def setUp(self):
        self.interactivetTest = False
    # interactivetTest = True

    @attr("LAR")
    @unittest.skipIf(not ("larcc" in sys.modules), "larcc is not installed")
    def test_vessel_tree_lar(self):
        import teigen.tb_lar
        tvg = TreeBuilder(teigen.tb_lar.TBLar)
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        output = tvg.buildTree() # noqa
        if self.interactiveTests:
            tvg.show()

    def test_vessel_tree_vtk(self):
        tvg = TreeBuilder('vtk')
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        output = tvg.buildTree() # noqa
        tvg.show()
        tvg.saveToFile("tree_output.vtk")

    @unittest.skipIf(not ("skelet3d" in sys.modules), "skelet3d is not installed")
    def test_vessel_tree_vtk_from_skeleton(self):

        import skelet3d
        import skelet3d.skeleton_analyser
        import shutil

        fn_out = 'tree.vtk'
        if os.path.exists(fn_out):
            os.remove(fn_out)

        volume_data = np.zeros([3, 7, 9], dtype=np.int)
        volume_data [:, :, 1:3] = 1
        volume_data [:, 5, 2:9] = 1
        volume_data [:, 0:7, 5] = 1
        skelet = skelet3d.skelet3d(volume_data)

        skan = skelet3d.skeleton_analyser.SkeletonAnalyser(skelet, volume_data=volume_data, voxelsize_mm=[1,1,1])
        stats = skan.skeleton_analysis()

        tvg = TreeBuilder('vtk')
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        tvg.tree_data = stats
        output = tvg.buildTree() # noqa
        tvg.saveToFile(fn_out)
        os.path.exists(fn_out)

    # TODO finish this test
    def test_vessel_tree_vol(self):
        import teigen.tb_volume
        tvg = TreeBuilder(teigen.tb_volume.TBVolume)
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        output = tvg.buildTree() # noqa
        # tvg.show()
        # if self.interactiveTests:
        #     tvg.show()

    def test_import_new_vt_format(self):

        tvg = TreeBuilder()
        yaml_path = os.path.join(path_to_script, "vt_biodur.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [150, 150, 150]
        data3d = tvg.buildTree()

    def test_cylinders_generator(self):
        from teigen.generators.cylinders import CylinderGenerator

        cg = CylinderGenerator()
        cg.run()


    def test_get_line_nodes(self):
        import teigen.geometry3d as g3
        nodes = g3.get_points_in_line_segment([10, 13, 22], [1, 13, 22], 3)
        expected_x = [10, 7, 4, 1]

        self.assertAlmostEqual(nodes[1][0], expected_x[1])
        self.assertAlmostEqual(nodes[2][0], expected_x[2])
        self.assertAlmostEqual(nodes[3][0], expected_x[3])

    def test_tree_generator(self):
        import numpy as np
        tree_data = {

        }
        element_number = 10
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
                    "nodeA_ZYX_mm": pts[i-1],
                    "nodeB_ZYX_mm": pts[i],
                    #         "nodeB_ZYX_mm": np.random.random(3) * 100,
                    "radius_mm": 1
                }
                tree_data[i+length] = edge

        tvg = TreeBuilder('vtk')
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        # tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        tvg.tree_data = tree_data
        output = tvg.buildTree() # noqa
        # tvg.show()
        tvg.saveToFile("test_tree_output.vtk")


        tvgvol = TreeBuilder('vol')
        tvgvol.voxelsize_mm = [1, 1, 1]
        tvgvol.shape = [100, 100, 100]
        tvgvol.tree_data = tree_data
        outputvol = tvgvol.buildTree()
        tvgvol.saveToFile("tree_volume.pklz")
        # self.assertTrue(False)



    def test_io3d(self):
        import io3d
        data3d = np.zeros([10,10,10])
        segmentation = np.zeros([10,10,10])

        data3d [2:7,:3:5, :6] = 100
        datap = {
            "data3d": data3d,
            # "segmentation": segmentation,
            "voxelsize_mm": [1,1,1]
        }
        io3d.write(datap, "file1.pklz")

    def test_skimage_io_imsave(self):
        import skimage.io
        data3d = np.zeros([10,10,10])
        segmentation = np.zeros([10,10,10])

        data3d [2:7,:3:5, :6] = 100
        datap = {
            "data3d": data3d,
            # "segmentation": segmentation,
            "voxelsize_mm": [1,1,1]
        }
        skimage.io.imsave("skiamge.png", data3d[0])


if __name__ == "__main__":
    unittest.main()
