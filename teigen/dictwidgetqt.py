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
try:
    QString = unicode
except NameError:
    # Python 3
    QString = str

logger = logging.getLogger(__name__)

import argparse


import PyQt5


from PyQt5.QtWidgets import (QLabel, QPushButton, QApplication, QGridLayout,
                         QSpinBox, QLineEdit, QCheckBox, QDoubleSpinBox,
                         QBoxLayout, QRadioButton, QWidget, QComboBox,
                         QBoxLayout)


from PyQt5 import QtGui, QtWidgets
import sys
import numpy as np

from pyqtconfig import ConfigManager
import copy

# class ConfigManager():
#     # TODO check config manager
#     def __init__(self):
#         self.config = {}
#         self.defaults = {}
#         pass
#
#     def set_defaults(self, defaults):
#         self.defaults = defaults
#         self.config.update(defaults)
#
#     def add_handler(self, key, widget):
#         import pdb; pdb.set_trace()
#         value = widget.value()
#         self.config[key] = value


class DictWidget(QtWidgets.QWidget):

    def __init__(self, config_in, ncols=2, captions=None, hide_keys=None,
                 horizontal=False, show_captions=True, accept_button=False,
                 config_manager=None, radiobuttons=None, dropdownboxes=None):
        """
        :param config_in:  dictionary
        :param ncols:
        :param captions:
        :param radiobuttons: {"name_of_radiobutton": [
                                 ["orig_caption_1", orig_caption_2],
                                 default_value]
                                 ]
                            }
        """
        super(DictWidget, self).__init__()
        if captions is None:
            captions = {}
        if hide_keys is None:
            hide_keys = []
        self.config_in = config_in
        self.ncols = ncols
        self.captions = captions
        self.accept_button = accept_button
        self.hide_keys = copy.copy(hide_keys)
        self.horizontal = horizontal
        self.show_captions = show_captions
        if radiobuttons is None:
            radiobuttons = {}
        self.radiobuttons = radiobuttons
        if dropdownboxes is None:
            dropdownboxes = {}
        self.dropdownboxes = dropdownboxes

        # hide also temp keys for lists and ndarrays
        # due to load default params
        self._get_tmp_composed_keys(config_in)
        self.hide_keys.extend(self._tmp_composed_keys_list)

        if config_manager is None:
            self.config = ConfigManager()
            self.config.set_defaults(config_in)
        else:
            self.config = config_manager
        self.mainLayout = QGridLayout(self)
        self.setLayout(self.mainLayout)
        self.widgets = {}
        self.grid_i = 0
        self.init_ui()

    def _get_tmp_composed_keys(self, cfg):
        """
        Get list of temporary keys for lists and ndarrays. This keys are used for reconstruction.

        vytvoří to seznam pomocných klíčů pro seznamy a ndarray
        :param cfg:
        :return:
        """
        self._tmp_composed_keys_dict = {}
        self._tmp_composed_keys_list = []
        toappend = {}
        for key, value in cfg.items():
            if key in self.dropdownboxes.keys():
                continue
            if key in self.radiobuttons.keys():
                continue
            if type(value) in (list, np.ndarray):
                self._tmp_composed_keys_dict[key] = []
                array = np.asarray(value)
                key_array_i = 0
                for val in array.tolist():
                    # key_i = (key, key_array_i)
                    key_i = ComposedDictMetadata((key, key_array_i))
                    self._tmp_composed_keys_dict[key].append(key_i)
                    self._tmp_composed_keys_list.append(key_i)
                    key_array_i += 1
                    toappend[key_i] = val
        cfg.update(toappend)

    def _create_dropdownbox(self, key, value):
        # atomic_widget = QWidget()
        row, col = self.__calculate_new_grid_position()
        # layout = QBoxLayout(self.horizontal)
        atomic_widget = QComboBox()
        # atomic_widget.addItem("C")
        # atomic_widget.addItem("C++")
        values = self.dropdownboxes[key]
        # values = self.dropdownboxes[key][0]
        if value is not None and value in values:
            # vali = atomic_widget.findText(value)
            atomic_widget.findText(value)
        atomic_widget.addItems(values)
        # this does not work. I used findText()
        # atomic_widget.setCurrentIndex(vali)
        # layout.addWidget(cb)

        # atomic_widget.setLayout(layout)
        return atomic_widget

    def _create_radiobutton(self, key, value):
        atomic_widget = QWidget()
        layout = QBoxLayout(self.horizontal)
        for i, rbkey in enumerate(self.radiobuttons[key][0]):
            b1 = QRadioButton("Button1")
            if i == self.radiobuttons[key][1]:
                b1.setChecked(True)
            # b1.toggled.connect(lambda:self.btnstate(self.b1))
            layout.addWidget(b1)
        atomic_widget.setLayout(layout)
        return atomic_widget

    def init_ui(self):
        # self.widgets = {}
        # self.grid_i = 0
        grid = self.mainLayout

        for key, value in self.config_in.items():

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
            btn_accept.clicked.connect(self.btn_accept)
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

        if key in self.dropdownboxes.keys():
            atomic_widget = self._create_dropdownbox(key,value)
            self.config.add_handler(key, atomic_widget)
        elif key in self.radiobuttons.keys():
            atomic_widget = self._create_radiobutton(key,value)
            self.config.add_handler(key, atomic_widget)
        elif type(value) is int:
            atomic_widget = QSpinBox()
            atomic_widget.setRange(-100000, 100000)
            self.config.add_handler(key, atomic_widget)
        elif type(value) is float:
            atomic_widget = QDoubleSpinBox()
            atomic_widget.setDecimals(6)
            atomic_widget.setMaximum(1000000000)
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

    def __calculate_new_grid_position(self):
        row = self.grid_i / self.ncols
        col = (self.grid_i % self.ncols) * 2
        self.grid_i += 1
        if self.horizontal:
            return col, row
        return row, col

    def btn_accept(self):
        logger.debug(("btn_accept: " + str(self.config_as_dict())))

    def on_config_update(self):
        pass

    def config_as_dict(self):
        def _primitive_type(value):
            if type(value) == QString:
                value = str(value)
            return value

        dictionary = self.config.as_dict()
        dictionary = copy.copy(dictionary)
        for key, value in dictionary.items():
            from PyQt5.QtCore import pyqtRemoveInputHook
            #pyqtRemoveInputHook()
            # import ipdb; ipdb.set_trace() #  noqa BREAKPOINT
            # if type(key) == tuple:
            if type(key) == ComposedDictMetadata:
                dict_key, value_index = key
                dict_key = (dict_key)
                # if dict_key not in dictionary.keys():
                #     dictionary[dict_key] = {}
                dictionary[dict_key][value_index] = _primitive_type(value)

            else:
                dictionary[key] = _primitive_type(value)

        for key in dictionary.keys():
            # if type(key) == tuple:
            if type(key) == ComposedDictMetadata:
                dictionary.pop(key)

        # for key, value in dictionary.items():
        #     if type(key) == tuple:
        #         dictionary

        return dictionary


def complicated_to_yaml(cfg):
    """
    write complex dict structure to yaml
    :param cfg:
    :return:
    """
    from ruamel.yaml import YAML
    yaml = YAML()
    # convert values to json
    isconverted = {}
    for key, value in cfg.items():
        if type(value) in (str, int, float, bool):

            isconverted[key] = False
            if type(value) is str:
                pass

        else:
            isconverted[key] = True
            cfg[key] = yaml.dump(value, default_flow_style=True)
    return cfg


class ComposedDictMetadata(tuple):
    pass


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
    cfg = {"bool": True, "int": 5, 'str': 'strdrr'}
    captions = {"int": "toto je int"}
    cw = DictWidget(cfg, captions=captions)
    cw.show()
    app.exec_()


if __name__ == "__main__":
    main()
