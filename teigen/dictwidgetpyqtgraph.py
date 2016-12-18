#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © %YEAR%  <>
#
# Distributed under terms of the %LICENSE% license.

"""

"""

import logging

logger = logging.getLogger(__name__)
import argparse

import PyQt4
from PyQt4.QtGui import QGridLayout, QLabel,\
    QPushButton, QLineEdit, QApplication, QWidget, QGridLayout, QSpinBox, QLineEdit, QCheckBox,\
        QComboBox, QTextEdit, QDialog, QMainWindow, QDoubleSpinBox

from PyQt4 import QtGui
import sys
import os.path
import copy
import numpy as np

from pyqtconfig import ConfigManager

import inspect
import collections
import numpy as np

import pyqtgraph.parametertree.parameterTypes as pTypes

def to_pyqtgraph_struct(name, value, opts={}):
    """
    Prepare structure for visualization by pyqtgraph tree.

    :param name:
    :param value:
    :param opts:
    :return:
    """

    tp = value.__class__.__name__
    if tp in (
        'list', 'ndarray', 'OrderedDict', 'dict',
        'int', 'float', 'bool', 'str', 'color', 'colormap'
    ):
        pass
    else:
        # some custome object
        return value
    item_properties = {
        "name": name,
        'value': value,
        'type': tp,
    }
    # if key in params.keys():

    children_properties = {}
    if "children" in opts.keys():
        children_properties = opts.pop('children')

    item_properties.update(opts)
    item_properties['reconstruction_type'] = tp

    if tp in ('list', 'ndarray', 'OrderedDict', 'dict'):
        # key_parameters['type'] = key_parameters['type']
        item_properties['type'] = 'group'
        item_properties.pop('value')
        if tp == 'list':
            children_key_value = collections.OrderedDict(zip(map(str, range(len(value))), value))
        elif (tp == 'ndarray'):
            value_list = value.to_list()
            children_key_value = collections.OrderedDict(zip(map(str, range(len(value_list))), value_list))
        elif tp in ('dict', 'OrderedDict'):
            children_key_value = value

        children_list = []
        for keyi, vali in children_key_value.items():
            children_properties_i = {}
            if keyi in children_properties:
                children_properties_i.update(children_properties[keyi])
            children_item = to_pyqtgraph_struct(keyi, vali, children_properties_i)
            children_list.append(children_item)

        item_properties['children'] = children_list

            # value = value_list
            # print key_parameters
    return item_properties

    #     if ntype is not None:
    #         key_parameters['type'] = ntype
    #         outdict.append(key_parameters)
    #
    # return outdict

def from_pyqtgraph_struct(dct):
    output = {}
    key = dct['name']
    if 'children' in dct.keys():
        reconstruction_type = 'dict'

        if 'reconstruction_type' in dct.keys():
            reconstruction_type = dct['reconstruction_type']

        if reconstruction_type in ('list', 'ndarray'):
            children_list = []
            for child in dct['children']:
                keyi, valuei = from_pyqtgraph_struct(dct['children'][child])
                children_list.append(valuei)
            value = children_list

        elif reconstruction_type in ('dict', 'OrderedDict'):
            children_dict = {}
            for child in dct['children']:
                child_item = dct['children'][child]
                keyi, valuei = from_pyqtgraph_struct(child_item)
                children_dict[keyi] = valuei
            value = children_dict


    else:
        value = dct['value']

    return key, value



class ListParameter(pTypes.GroupParameter):
    """
    New keywords

    titles: is list of titles
    value: list of values
    """
    def __init__(self, **opts):
        values = opts.pop('value')
        parent_opts={
            'name': opts.pop('name'),
            'type': 'bool',
            'value': values
        }
        if 'title' in opts.keys():
            parent_opts['title']  = opts.pop('title')
        if "reconstruction_type" in opts.keys():
            parent_opts['reconstruction_type']  = opts.pop('reconstruction_type')
        # opts['type'] = 'bool'
        # opts['value'] = True
        pTypes.GroupParameter.__init__(self, **parent_opts)
        # gp = pTypes.GroupParameter(name=opts['name'], title=opts['title'])
        if 'names' in opts.keys():
            names = opts['names']
        else:
            names =  map(str, range(len(values)))

        for i in range(len(values)):
            opts['name'] = names[i]
            opts['value'] = values[i]
            self.addChild(opts)

        for child in self.childs:
            child.sigValueChanged.connect(self.valuesChanged)

        self.sigValueChanged.connect(self.valueChanged)

    def valuesChanged(self):
        new_val = []
        for i in range(len(self.childs)):
            new_val.append(self.childs[i].value())

        self.setValue(new_val)
        print new_val

    def valueChanged(self):
        new_val = self.value()
        for i in range(len(self.childs)):
            self.childs[i].setValue(new_val[i])
        print new_val



class AreaSamplingParameter(pTypes.GroupParameter):
    def __init__(self, **opts):
        opts['type'] = 'bool'
        opts['value'] = True
        voxelsize_mm = [1.0, 1.0, 1.0]
        areasize_mm = [10.0, 10.0, 10.0]
        areasize_px = [10, 10, 10]

        if "voxelsize_mm" in opts.keys():
            voxelsize_mm = opts.pop('voxelsize_mm')
        if "areasize_mm" in opts.keys():
            areasize_mm = opts.pop('areasize_mm')
        if "areasize_px" in opts.keys():
            areasize_px = opts.pop('areasize_px')

        pTypes.GroupParameter.__init__(self, **opts)

        self.p_voxelsize_mm = ListParameter(name="voxelsize_mm", value=voxelsize_mm, type='float', suffix='mm', siPrefix=False, reconstruction_type='list')
        self.p_areasize_px = ListParameter(name="areasize_px", value=areasize_px, type='int', suffix='px', siPrefix=False, reconstruction_type='list')
        self.p_areasize_mm = ListParameter(name="areasize_mm", value=areasize_mm, type='float', suffix='mm', siPrefix=False, reconstruction_type='list')

        self.addChild(self.p_voxelsize_mm)
        self.addChild(self.p_areasize_mm)
        self.addChild(self.p_areasize_px)

        self.p_areasize_mm.sigValueChanged.connect(self.areasize_mChanged)
        self.p_areasize_px.sigValueChanged.connect(self.areasize_pxChanged)


    def areasize_mChanged(self):
        as_m = np.asarray(self.p_areasize_mm.value())
        vs_m = np.asarray(self.p_voxelsize_mm.value()).astype(np.float)
        val = (as_m / vs_m).astype(np.int).tolist()
        self.p_areasize_px.setValue(
            val,
            blockSignal=self.areasize_pxChanged)

    def areasize_pxChanged(self):
        val = (np.asarray(self.p_voxelsize_mm.value()) * np.asarray(self.p_areasize_px.value())).tolist()
        self.p_areasize_mm.setValue(
            val,
            blockSignal=self.areasize_mChanged)

    def voxelsizeChanged(self):
        self.z_size_px.setValue(int(self.z_size_m.value() / self.z_m.value()), blockSignal=self.z_size_pxChanged)
    # def z_mChanged(self):
    #     self.z_.setValue(1.0 / self.a.value(), blockSignal=self.bChanged)
    def z_size_mChanged(self):
        self.z_size_px.setValue(int(self.z_size_m.value() / self.z_m.value()), blockSignal=self.z_size_pxChanged)
    def z_size_pxChanged(self):
        self.z_size_m.setValue(int(self.z_size_px.value() * self.z_m.value()), blockSignal=self.z_size_mChanged)

    # def aChanged(self):
    #     self.b.setValue(1.0 / self.a.value(), blockSignal=self.bChanged)
    #
    # def bChanged(self):
    #     self.a.setValue(1.0 / self.b.value(), blockSignal=self.aChanged)

def main():
    logger = logging.getLogger()

    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    logger.addHandler(ch)

    # create file handler which logs even debug messages
    # fh = logging.FileHandler('log.txt')
    # fh.setLevel(logging.DEBUG)
    # formatter = logging.Formatter(
    #     '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # fh.setFormatter(formatter)
    # logger.addHandler(fh)
    # logger.debug('start')

    # input parser
    parser = argparse.ArgumentParser(
        description=__doc__
    )
    parser.add_argument(
        '-i', '--inputfile',
        default=None,
        required=True,
        help='input file'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='Debug mode')
    args = parser.parse_args()

    if args.debug:
        ch.setLevel(logging.DEBUG)


    app = QApplication(sys.argv)
    cfg = {"bool": True, "int":5, 'str': 'strdrr'}
    captions = {"int": "toto je int"}
    to_pyqtgraph_struct(cfg)
    # cw = DictWidget(cfg, captions=captions)
    # cw.show()
    # app.exec_()


if __name__ == "__main__":
    main()
