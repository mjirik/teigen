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

from dili import get_default_args, subdict

class DictWidget(QtGui.QWidget):
    def __init__(self, config_in, ncols=2, captions={}, hide_keys=[], horizontal=False, show_captions=True, accept_button=False, config_manager=None):
        """

        :param config_in:  dictionary
        :param ncols:
        :param captions:
        """
        super(DictWidget, self).__init__()
        self.config_in = config_in
        self.ncols = ncols
        self.captions = captions
        self.accept_button = accept_button
        self.hide_keys = copy.copy(hide_keys)
        self.horizontal = horizontal
        self.show_captions = show_captions

        # hide also temp keys for lists and ndarrays
        # due to load default params
        self._get_tmp_composed_keys(config_in)
        rr = self.hide_keys.extend(self._tmp_composed_keys_list)

        if config_manager is None:
            self.config = ConfigManager()
            self.config.set_defaults(config_in)
        else:
            self.config = config_manager
        self.init_ui()

    def _get_tmp_composed_keys(self, cfg):
        # vytvoří to seznam pomocných klíčů pro seznamy a ndarray
        self._tmp_composed_keys_dict = {}
        self._tmp_composed_keys_list = []
        for key, value in cfg.iteritems():
            if type(value) in (list, np.ndarray):
                self._tmp_composed_keys_dict[key] = []
                array = np.asarray(value)
                key_array_i = 0
                for val in array.tolist():
                    # key_i = key + str(hgrid_i)
                    key_i = (key, key_array_i)
                    self._tmp_composed_keys_dict[key].append(key_i)
                    self._tmp_composed_keys_list.append(key_i)
                    key_array_i += 1
                    cfg[key_i] = val

        # self._tmp_composed_keys.keys()




    def complicated_to_yaml(self, cfg):
        import yaml
        # convert values to json
        isconverted = {}
        for key, value in cfg.iteritems():
            if type(value) in (str, int, float, bool):

                isconverted[key] = False
                if type(value) is str:
                    pass

            else:
                isconverted[key] = True
                cfg[key] = yaml.dump(value, default_flow_style=True)
        return cfg


    def init_ui(self):
        self.mainLayout = QGridLayout(self)
        self.widgets = {}
        grid = self.mainLayout
        self.grid_i = 0


        for key, value in self.config_in.iteritems():

            if key in self.hide_keys:
                continue
            if key in self.captions.keys():
                caption = self.captions[key]
            else:
                caption = key

            atomic_widget = self.__get_widget_for_primitive_types(key, value)
            if atomic_widget is None:
                if type(value) in (list, np.ndarray):
                    array = np.asarray(value)
                    atomic_widget = self._create_sub_grid_from_ndarray(key, array)
                    # dc = dict(zip(range(len(vl)), list(vl.astype(str))))
                    # atomic_widget = DictWidget(config_in=dc, show_captions=False, horizontal=True, config_manager=self.config)
                    # atomic_widget.show()
                    row, col = self.__calculate_new_grid_position()
                    grid.addWidget(QLabel(caption), row, col + 1)
                    grid.addLayout(atomic_widget, row, col + 2)
                    continue
                else:
                    logger.error("Unexpected type in config dictionary")

                continue

            # import ipdb; ipdb.set_trace()
            row, col = self.__calculate_new_grid_position()
            grid.addWidget(QLabel(caption), row, col + 1)
            grid.addWidget(atomic_widget, row, col + 2)

        # gd.setColumnMinimumWidth(text_col, 500)

        if self.accept_button:
            btn_accept = QPushButton("Accept", self)
            btn_accept.clicked.connect(self.btnAccept)
            text_col = (self.ncols * 2) + 3
            grid.addWidget(btn_accept, (self.grid_i / 2), text_col)

        self.config.updated.connect(self.on_config_update)

    # def __add_line
    def __get_widget_for_primitive_types(self, key, value):

        """
        return right widget and connect the value with config_manager
        :param key:
        :param value:
        :return:
        """

        if type(value) is int:
            atomic_widget = QSpinBox()
            atomic_widget.setRange(-100000, 100000)
            self.config.add_handler(key, atomic_widget)
        elif type(value) is float:
            atomic_widget = QDoubleSpinBox()
            self.config.add_handler(key, atomic_widget)
        elif type(value) is str:
            atomic_widget = QLineEdit()
            self.config.add_handler(key, atomic_widget)
        elif type(value) is bool:
            atomic_widget = QCheckBox()
            self.config.add_handler(key, atomic_widget)
        else:
            return None
        return atomic_widget

    def _create_sub_grid_from_ndarray(self, key, ndarray):
        hgrid = QGridLayout(self)
        hgrid_i = 0
        for val in ndarray.tolist():
            # key_i = key + str(hgrid_i)
            key_i = (key, hgrid_i)
            atomic_widget = self.__get_widget_for_primitive_types(key_i, val)

            hgrid.addWidget(atomic_widget, 0, hgrid_i)
            hgrid_i += 1

        return hgrid


    def __calculate_new_grid_position(self): #, atomic_widget, caption, grid):
        row = self.grid_i / self.ncols
        col = (self.grid_i % self.ncols) * 2
        self.grid_i += 1
        if self.horizontal:
             return col, row
        return row, col


    def btnAccept(self):
        print self.config_as_dict()

    def on_config_update(self):
        pass


    def config_as_dict(self):
        def _primitive_type(value):
            if type(value) == PyQt4.QtCore.QString:
                value = str(value)
            return value

        dictionary = self.config.as_dict()
        dictionary = copy.copy(dictionary)
        for key, value in dictionary.iteritems():
            from PyQt4.QtCore import pyqtRemoveInputHook
            pyqtRemoveInputHook()
            # import ipdb; ipdb.set_trace() #  noqa BREAKPOINT
            if type(key) == tuple:

                dict_key, value_index = key
                dict_key = (dict_key)
                # if dict_key not in dictionary.keys():
                #     dictionary[dict_key] = {}

                dictionary[dict_key][value_index] = _primitive_type(value)

            else:
                dictionary[key] = _primitive_type(value)

        for key in dictionary.keys():
            if type(key) == tuple:
                dictionary.pop(key)


        # for key, value in dictionary.iteritems():
        #     if type(key) == tuple:
        #         dictionary

        return dictionary


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
    cw = DictWidget(cfg, captions=captions)
    cw.show()
    app.exec_()


if __name__ == "__main__":
    main()