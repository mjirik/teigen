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

def dict_to_pyqtgraph(key_value={}, manual_parameters={}):

    outdict = []

    for key in key_value:
        value = key_value[key]
        key_parameters = {
            "name":key,
            'value': value,
        }
        ntype = None
        tp = type(value)
        if tp == int:
            ntype = 'int'
        elif (tp == bool):
            ntype = 'bool'
        elif (tp == str):
            ntype = 'str'
        elif (tp == float):
            ntype = 'float'
        elif tp in (list, np.ndarray, collections.OrderedDict):
            ntype = 'group'
            key_parameters.pop('value')
            if tp == list:
                children_key_value = collections.OrderedDict(zip(map(str, range(len(value))), value))
            elif (tp == np.ndarray):
                value_list = value.to_list()
                children_key_value = collections.OrderedDict(zip(map(str, range(len(value_list))), value_list))
            elif tp in (dict, collections.OrderedDict):
                children_key_value = value
            children_dict = dict_to_pyqtgraph(key_value=children_key_value)

            key_parameters['children'] = children_dict
            print key_parameters


        if ntype is not None:
            key_parameters['type'] = ntype
            outdict.append(key_parameters)

    return outdict

# def add_item()


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