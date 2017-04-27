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

    @attr('interactive')
    def test_qt_file_dialog(self):
        import PyQt4
        from PyQt4.QtGui import QApplication, QFileDialog
        app = QApplication(sys.argv)
        # fname = QFileDialog.getOpenFileName(None, 'Open file',
        #  'c:\\',"Image files (*.jpg *.gif)")
        filename = QFileDialog.getSaveFileName(
            None,
            "Save file",
            'c:\\ahojfile%06d.jpg',
            ""
            # "Image files (*.jpg *.gif)"
            # ""
        )

        # QFileDialog.getExistingDirectory()

    @attr('interactive')
    def test_qt_dictwidget(self):
        import PyQt4
        from PyQt4.QtGui import QApplication, QFileDialog
        from teigen.dictwidgetqt import DictWidget
        app = QApplication(sys.argv)
        cfg = {"bool": True, "int":5, 'str': 'strdrr', 'vs':[1.0, 2.5, 7]}
        captions = {"int": "toto je int"}
        cw = DictWidget(cfg, captions=captions)
        cw.show()
        app.exec_()

    @attr('interactive')
    def test_qt_cylinderwidget(self):
        import PyQt4
        from PyQt4.QtGui import QApplication, QFileDialog
        from teigen.dictwidgetqt import DictWidget
        from teigen.generators.cylindersqt import CylindersWidget
        app = QApplication(sys.argv)
        cw = CylindersWidget()
        cw.show()
        app.exec_()

    # tenhle by měl fungovat i neinteraktivně
    # @attr('interactive')
    @unittest.skip('some VTK problem with CylinderWidget.run()')
    def test_qt_cylinderwidget_run(self):
        import PyQt4
        from PyQt4.QtGui import QApplication, QFileDialog
        from teigen.dictwidgetqt import DictWidget
        from teigen.generators.cylindersqt import CylindersWidget
        app = QApplication(sys.argv)
        cw = CylindersWidget()
        cw.show()
        cw.run()
        cw.deleteLater()
        app.deleteLater()
        # app.exec_()

    @attr('interactive')
    def test_qt_set_dir_widget_interactive(self):
        import PyQt4
        from PyQt4.QtGui import QApplication, QFileDialog
        from teigen.dictwidgetqt import DictWidget
        from teigen.iowidgetqt import SetDirWidget
        app = QApplication(sys.argv)
        cfg = {"bool": True, "int":5, 'str': 'strdrr', 'vs':[1.0, 2.5, 7]}
        captions = {"int": "toto je int"}
        cw = SetDirWidget("~/lisa_data", "output dir")
        cw.show()
        app.exec_()


    # @unittest.skip('some VTK problem')
    @attr('interactive')
    def test_qt_set_dir_widget(self):
        import PyQt4
        from PyQt4.QtGui import QApplication, QFileDialog
        from teigen.dictwidgetqt import DictWidget
        from teigen.iowidgetqt import SetDirWidget
        app = QApplication(sys.argv)
        cw = SetDirWidget("~/lisa_data", "output dir")
        cw.show()
        cw.deleteLater()
        app.deleteLater()


def dist_to_vectors(v1, vlist):
    import numpy as np
    out = []
    for v2 in vlist:
        dist = np.linalg.norm(v1-v2)
        out.append(dist)
    return np.asarray(out)




if __name__ == "__main__":
    unittest.main()
