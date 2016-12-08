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

def dict_to_pyqtgraph(key, value, params={}):

    tp = value.__class__.__name__
    key_parameters = {
        "name":key,
        'value': value,
        'type': tp,
    }
    # if key in params.keys():

    children_params = {}
    if "children" in params.keys():
        children_params = params.pop('children_params')

    key_parameters.update(params)
    key_parameters['reconstruction_type'] = tp

    # ntype = None
    # if tp == int:
    #     ntype = 'int'
    # elif (tp == bool):
    #     ntype = 'bool'
    # elif (tp == str):
    #     ntype = 'str'
    # elif (tp == float):
    #     ntype = 'float'
    if tp in ('list', 'ndarray', 'OrderedDict'):
        # key_parameters['type'] = key_parameters['type']
        key_parameters['type'] = 'group'
        key_parameters.pop('value')
        if tp == 'list':
            children_key_value = collections.OrderedDict(zip(map(str, range(len(value))), value))
        elif (tp == 'ndarray'):
            value_list = value.to_list()
            children_key_value = collections.OrderedDict(zip(map(str, range(len(value_list))), value_list))
        elif tp in ('dict', 'OrderedDict'):
            children_key_value = value

        children_list = []
        for keyi, vali in children_key_value.items():
            children_item = dict_to_pyqtgraph(keyi, vali, children_params)
            children_list.append(children_item)

        key_parameters['children'] = children_list

            # value = value_list
            # print key_parameters
    return key_parameters

    #     if ntype is not None:
    #         key_parameters['type'] = ntype
    #         outdict.append(key_parameters)
    #
    # return outdict

# def add_item()
def pyqtgraph_to_dict(dct):
    output = {}
    key = dct['name']
    if 'children' in dct.keys():
        reconstruction_type = 'dict'
        # if reconsturuction type is list
        if reconstruction_type == 'dict':
            children_dict = {}
            for child in dct['children']:
                child_item = dct['children'][child]
                keyi, valuei = pyqtgraph_to_dict(child_item)
                children_dict[keyi] = valuei
            value = children_dict

        elif reconstruction_type == 'list':
            children_list = []
            for child in dct['children']:
                keyi, valuei = pyqtgraph_to_dict(dct['children'][child])
                children_list.append(valuei)
            value = children_list

    else:
        value = dct['value']

    return key, value




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
    dict_to_pyqtgraph(cfg)
    # cw = DictWidget(cfg, captions=captions)
    # cw.show()
    # app.exec_()


if __name__ == "__main__":
    main()