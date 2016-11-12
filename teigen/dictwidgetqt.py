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

from pyqtconfig import ConfigManager



class DictWidget(QtGui.QWidget):
    def __init__(self, config_in, ncols=2):
        super(DictWidget, self).__init__()
        self.config_in = config_in
        self.ncols = ncols
        self.config = ConfigManager()
        self.init_ui()


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
        gd = self.mainLayout

        gd_max_i = 0
        for key, value in self.config_in.iteritems():
            if type(value) is int:
                sb = QSpinBox()
                sb.setRange(-100000, 100000)
            elif type(value) is float:
                sb = QDoubleSpinBox()
            elif type(value) is str:
                sb = QLineEdit()
            elif type(value) is bool:
                sb = QCheckBox()
            else:
                logger.error("Unexpected type in config dictionary")


            row = gd_max_i / self.ncols
            col = (gd_max_i % self.ncols) * 2

            gd.addWidget(QLabel(key),row, col +1)
            gd.addWidget(sb, row, col + 2)
            self.config.add_handler(key, sb)
            gd_max_i += 1

        text_col = (self.ncols * 2) + 3
        gd.setColumnMinimumWidth(text_col, 500)

        btn_accept = QPushButton("Accept", self)
        btn_accept.clicked.connect(self.btnAccept)
        gd.addWidget(btn_accept, (gd_max_i / 2), text_col)

        self.config.updated.connect(self.on_config_update)

    def btnAccept(self):
        print self.config_as_dict()

    def on_config_update(self):
        pass


    def config_as_dict(self):
        dictionary = self.config.as_dict()
        for key, value in dictionary.iteritems():
            from PyQt4.QtCore import pyqtRemoveInputHook
            pyqtRemoveInputHook()
            # import ipdb; ipdb.set_trace() #  noqa BREAKPOINT
            if type(value) == PyQt4.QtCore.QString:
                value = str(value)
            dictionary[key] = value

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
    cfg = {"bool": True, "int":1}
    cw = DictWidget(cfg)
    cw.show()
    app.exec_()


if __name__ == "__main__":
    main()