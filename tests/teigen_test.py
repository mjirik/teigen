# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)
import unittest
import sys

from nose.plugins.attrib import attr
import teigen
import io3d

class MyTestCase(unittest.TestCase):

    @attr('interactive')
    def test_teigen_gui_interactive(self):
        import os.path as op
        params = io3d.misc.obj_from_file(op.expanduser("~/teigen_data/038/slice_parameters.yaml"))
        import PyQt4
        from PyQt4.QtGui import QApplication, QFileDialog
        # from teigen.dictwidgetqt import DictWidget
        import teigen
        import teigen.geometry3d
        import teigen.gui
        app = QApplication(sys.argv)
        cw = teigen.gui.TeigenWidget(config=params)
        cw.show()
        app.exec_()

    def test_teigen_gui(self):
        import PyQt4
        from PyQt4.QtGui import QApplication, QFileDialog
        # from teigen.dictwidgetqt import DictWidget
        import teigen
        import teigen.geometry3d
        import teigen.gui
        app = QApplication(sys.argv)
        cw = teigen.gui.TeigenWidget()
        cw.show()
        cw.deleteLater()
        app.deleteLater()

    def test_teigen(self):
        import teigen.gui
        tg = teigen.gui.Teigen()
        conf = {
            "areasampling":{
                "voxelsize_mm": [1., 1., 1.],
                "areasize_px": [150, 160, 170],
                "areasize_mm": [150, 160, 170],
            }
        }
        tg.update_config(**conf)
        tg.run()
        tg.save_volume()

if __name__ == '__main__':
    unittest.main()
