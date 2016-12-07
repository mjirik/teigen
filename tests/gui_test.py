#! /usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from nose.plugins.attrib import attr
import collections
import sys
import pyqtgraph
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType


class GuiTest(unittest.TestCase):
    @attr('interactive')
    def test_something(self):
        self.assertEqual(True, False)

    def test_pyqtgraph_import_dict(self):
        cfg = collections.OrderedDict(
            [["bool", True],
             ["int", 5],
             ['str', 'strdrr'],
             ['vs',[1.0, 2.5, 7]]]
        )
        captions = {"int": "toto je int"}
        import teigen.dictwidgetpyqtgraph
        params = teigen.dictwidgetpyqtgraph.dict_to_pyqtgraph(cfg)
        print params
        self.assertDictEqual(
            params[0],
            {'type': 'bool', 'name':'bool', 'value': True}
        )


    def test_pyqtgraph(self):
        cfg = collections.OrderedDict({"bool": True, "int":5, 'str': 'strdrr', 'vs':[1.0, 2.5, 7]})
        captions = {"int": "toto je int"}
        import teigen.dictwidgetpyqtgraph
        params = teigen.dictwidgetpyqtgraph.dict_to_pyqtgraph(cfg)
        print params

        from PyQt4.QtGui import QApplication, QFileDialog
        app = QApplication(sys.argv)
        p = Parameter.create(name='params', type='group', children=params)
        t = ParameterTree()
        print p.getValues()
        t.setParameters(p, showTop=False)
        t.show()

        app.exec_()

        self.assertTrue(False)


if __name__ == '__main__':
    unittest.main()
