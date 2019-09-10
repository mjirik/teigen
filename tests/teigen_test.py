#! /usr/bin/env python
# -*- coding: utf-8 -*-
import logging
logger = logging.getLogger(__name__)
import unittest
import sys

import pytest
# import teigen
# import io3d

import os.path as op
path_to_script = op.dirname(op.abspath(__file__))

class MyTestCase(unittest.TestCase):
    #
    @pytest.mark.interactive
    def test_teigen_gui_interactive(self):
        import os.path as op

        params = None
        # params = io3d.misc.obj_from_file(op.expanduser("~/teigen_data/038/slice_parameters.yaml"))

        from PyQt5.QtWidgets import QApplication, QFileDialog
        # from teigen.dictwidgetqt import DictWidget
        import teigen.gui

        app = QApplication(sys.argv)
        cw = teigen.gui.TeigenWidget(config=params)
        cw.show()
        app.exec_()

    @pytest.mark.interactive
    def test_teigen_gui_interactive_with_parameters(self):
        """
        reproduces undetected colision bug
        :return:
        """
        import os.path as op

        params = None
        # params = io3d.misc.obj_from_file(op.expanduser("~/teigen_data/038/slice_parameters.yaml"))
        params = {
            "generator_id": 3,
            "areasampling": {
                "voxelsize_mm": [1., 1., 1.],
                "areasize_px": [20, 20, 20],
                "areasize_mm": [20, 20, 20],
            },
            "postprocessing": {
                "measurement_resolution": 15,
                "measurement_multiplier": -1,
                "add_noise": False
            },
            "generators": {
                "Unconnected tubes": {
                    "element_number": 3,
                    "random_generator_seed": 110,
                    "radius_distribution_mean": 15,
                    "radius_distribution_maximum": 20,
                    "orientation_anisotropic": False,
                }
            }
        }
        # tg.update_config(**conf)


        from PyQt5.QtWidgets import QApplication
        # from teigen.dictwidgetqt import DictWidget
        import teigen.gui

        app = QApplication(sys.argv)
        cw = teigen.gui.TeigenWidget(use_default_config=True, config=params)
        cw.show()
        app.exec_()
    # def test_teigen_gui(self):
    #     import PyQt4
    #     from PyQt4.QtGui import QApplication, QFileDialog
    #     # from teigen.dictwidgetqt import DictWidget
    #     import teigen
    #     import teigen.geometry3d
    #     import teigen.gui
    #     app = QApplication(sys.argv)
    #     cw = teigen.gui.TeigenWidget()
    #     cw.show()
    #     cw.deleteLater()
    #     app.deleteLater()

    @pytest.mark.interactive
    def test_teigen_without_save(self):
        import teigen.gui

        tg = teigen.gui.Teigen()
        conf = {
            "generator_id": 3,
            "areasampling": {
                "voxelsize_mm": [1., 1., 1.],
                "areasize_px": [110, 120, 130],
                "areasize_mm": [110, 120, 130],
            },
            "postprocessing": {
                "measurement_multiplier": -1,
                "add_noise": False
            },
            "generators": {
                "Unconnected cylinders": {
                    "element_number": 10
                }
            }
        }
        tg.update_config(**conf)
        tg.step1()

    @pytest.mark.interactive
    def test_teigen_big(self):
        import teigen.gui

        tg = teigen.gui.Teigen()
        conf = {
            "areasampling": {
                "voxelsize_mm": [1., 1., 1.],
                "areasize_px": [210, 720, 730],
                "areasize_mm": [210, 720, 730],
            },
            "postprocessing": {
                "measurement_multiplier": -1,
                "add_noise": False
            },
            "generators": {
                "Unconnected cylinders": {
                    "element_number": 10
                }
            }
        }
        tg.update_config(**conf)
        tg.step1()
        tg.step2()

        # def test_teigen_small(self):
        #     import teigen.gui
        #     tg = teigen.gui.Teigen()
        #     conf = {
        #         "areasampling":{
        #             "voxelsize_mm": [1., 1., 1.],
        #             "areasize_px": [110, 120, 130],
        #             "areasize_mm": [110, 120, 130],
        #         },
        #         "postprocessing":{
        #             "measurement_multiplier":-1,
        #         }
        #     }
        #     tg.update_config(**conf)
        #     tg.run()
        #     tg.save_volume()

    def test_teigen_prepare_parameters_and_measurement(self):
        """
        Check string like generator_id
        :return:
        """
        print("test prepare parameters and measurement")
        import teigen.gui

        tg = teigen.gui.Teigen()
        tg.use_default_config()
        conf = {
            "generator_id": "Unconnected tubes",
            "areasampling": {
                "voxelsize_mm": [1., 1., 1.],
                "areasize_px": [110, 120, 130],
                "areasize_mm": [110, 120, 130],
            },
            "postprocessing": {
                # "measurement_multiplier": -1,
                "add_noise": False
            },
            "generators": {
                "Unconnected tubes": {
                    "element_number": 1
                }
            }
        }
        tg.update_config(**conf)
        tg.step1()
        params = tg.get_config_and_measurement()
        tg.step2()
        print(params)


    def test_teigen_read_tube_skeleton_from_file(self):
        """
        Read tube skeleton from file
        :return:
        """
        print("test read tube skeleton from file")
        import teigen.gui

        tg = teigen.gui.Teigen()
        tg.use_default_config()
        conf = {
            "generator_id": "Unconnected tubes",
            "areasampling": {
                "voxelsize_mm": [1., 1., 1.],
                "areasize_px": [110, 120, 130],
                "areasize_mm": [110, 120, 130],
            },
            "postprocessing": {
                # "measurement_multiplier": -1,
                "add_noise": False
            },
            "generators": {
                "Unconnected tubes": {
                    "element_number": 1
                }
            }
        }
        tg.update_config(**conf)
        tg.step1_by_load_tube_skeleton(
            op.join(path_to_script, "data_vt.yaml" ))
            #op.join(path_to_script, "vt_biodur.yaml" ))
        params = tg.get_config_and_measurement()
        tg.step2()
        print(params)

if __name__ == '__main__':
    unittest.main()
