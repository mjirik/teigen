# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)
import unittest
import sys

from nose.plugins.attrib import attr
import teigen

class MyTestCase(unittest.TestCase):

    @attr('interactive')
    def test_teigen_interactive(self):
        import PyQt4
        from PyQt4.QtGui import QApplication, QFileDialog
        # from teigen.dictwidgetqt import DictWidget
        import teigen
        import teigen.geometry3d
        import teigen.gui
        app = QApplication(sys.argv)
        cw = teigen.gui.TeigenWidget()
        cw.show()
        app.exec_()

    def test_teigen(self):
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
if __name__ == '__main__':
    unittest.main()
