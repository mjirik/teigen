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

def to_pyqtgraph_struct(name, value, properties={}):
    """
    Prepare structure for visualization by pyqtgraph tree.

    :param name:
    :param value:
    :param properties:
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
    if "children" in properties.keys():
        children_properties = properties.pop('children')

    item_properties.update(properties)
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



class ComplexParameter(pTypes.GroupParameter):
    def __init__(self, **opts):
        opts['type'] = 'bool'
        opts['value'] = True
        pTypes.GroupParameter.__init__(self, **opts)

        # self.addChild({'name': 'A = 1/B', 'type': 'float', 'value': 7, 'suffix': 'Hz', 'siPrefix': True})
        # self.addChild({'name': 'B = 1/A', 'type': 'float', 'value': 1/7., 'suffix': 's', 'siPrefix': True})
        # self.a = self.param('A = 1/B')
        # self.b = self.param('B = 1/A')
        # self.a.sigValueChanged.connect(self.aChanged)
        # self.b.sigValueChanged.connect(self.bChanged)
        pvoxelsize = pTypes.GroupParameter(name="voxelsize_mm", title="voxelsize [mm]")
        pvsz = pTypes.Parameter(name='z_m', type='float', value=0.01, suffix='m', siPrefix=True)
        pvsx = pTypes.Parameter(name='x_m', type='float', value=0.01, suffix='m', siPrefix=True)
        pvsy = pTypes.Parameter(name='y_m', type='float', value=0.01, suffix='m', siPrefix=True)
        pvoxelsize.addChild(pvsz)
        pvoxelsize.addChild(pvsx)
        pvoxelsize.addChild(pvsy)

        pareasize = pTypes.GroupParameter(name="shape", title="area size [px]")
        pasz = pTypes.Parameter(name='z', type='float', value=0.01, suffix='px', siPrefix=False)
        pasx = pTypes.Parameter(name='x', type='float', value=0.01, suffix='px', siPrefix=False)
        pasy = pTypes.Parameter(name='y', type='float', value=0.01, suffix='px', siPrefix=False)
        pareasize.addChild(pasz)
        pareasize.addChild(pasx)
        pareasize.addChild(pasy)
        # pvoxelsize.addChild({'name': 'z_m', 'type': 'float', 'value': 0.01, 'suffix': 'm', 'siPrefix': True})

        self.addChild({'name': 'z_m', 'type': 'float', 'value': 0.01, 'suffix': 'm', 'siPrefix': True})
        self.addChild({'name': 'z_size_px', 'type': 'int', 'value': 100, 'suffix': 'px', 'siPrefix': False})
        self.addChild({'name': 'z_size_m', 'type': 'float', 'value': 10.0, 'suffix': 'm', 'siPrefix': True})
        self.addChild({'name': 'x_m', 'type': 'float', 'value': 1/7., 'suffix': 's', 'siPrefix': True})
        self.addChild({'name': 'y_m', 'type': 'float', 'value': 1/7., 'suffix': 's', 'siPrefix': True})


        self.z_m = self.param('z_m')
        self.z_size_m = self.param('z_size_m')
        self.z_size_px = self.param('z_size_px')

        self.z_size_m.sigValueChanged.connect(self.z_size_mChanged)
        self.z_size_px.sigValueChanged.connect(self.z_size_pxChanged)

        self.addChild(pvoxelsize)
        self.addChild(pareasize)

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
