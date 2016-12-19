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
        params = teigen.dictwidgetpyqtgraph.to_pyqtgraph_struct(cfg)
        print params
        self.assertDictEqual(
            params[0],
            {'type': 'bool', 'name':'bool', 'value': True}
        )


    def test_pyqtgraph_from_manual_dict(self):
        params =  [
            {'type': 'int', 'name': 'int', 'value': 5},
            {'type': 'group', 'name': 'vs', 'children': [
                {'type': 'int', 'name': "2", 'value': 7}
            ]},
            {'type': 'bool', 'name': 'bool', 'value': True},
            {'type': 'str', 'name': 'str', 'value': 'strdrr'}
        ]

        p = Parameter.create(name='params', type='group', children=params)

    @attr('interactive')
    def test_pyqtgraph(self):
        cfg = collections.OrderedDict({
            "bool": True,
            "int":5,
            'str': 'strdrr',
            'vs':[1.0, 2.5, 7]
            # 'Area Sampling' : dictwidgetpyqtgraph.AreaSamplingParameter(name='Area Sampling')
        })
        captions = {"int": "toto je int"}

        opts = {}
        # opts = {
        #     "children": {
        #         "voxelsize_mm": {
        #             "title": 'voxelsize [mm]',
        #             "children": {
        #                 "0": {
        #                     "title": "z",
        #                     'suffix': 'm',
        #                     'siPrefix': True
        #                 },
        #                 "1": {
        #                     "title": "x",
        #                     'suffix': 'm',
        #                     'siPrefix': True
        #                 },
        #                 "2": {
        #                     "title": "y",
        #                     'suffix': 'm',
        #                     'siPrefix': True
        #                 }
        #             }
        #         }
        #     }
        # }
        import teigen.dictwidgetpyqtgraph
        params = teigen.dictwidgetpyqtgraph.to_pyqtgraph_struct('params', cfg, opts=opts)
        params['children'].append(
            teigen.dictwidgetpyqtgraph.AreaSamplingParameter(name='Area Sampling'))
        # print params

        # params[0]['title'] = "Pokusny title"
        # params[0]['my note'] = "poznamka"


        from PyQt4.QtGui import QApplication, QFileDialog
        app = QApplication(sys.argv)
        p = Parameter.create(**params)
        # p = Parameter.create(name='params', type='group', children=params)
        t = ParameterTree()
        print p.getValues()
        lst = p.saveState()
        vals = p.getValues()

        name, dict_again = teigen.dictwidgetpyqtgraph.from_pyqtgraph_struct(lst)
        t.setParameters(p, showTop=False)
        t.show()

        app.exec_()

        # self.assertTrue(False)


if __name__ == '__main__':
    unittest.main()
