#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © %YEAR%  <>
#
# Distributed under terms of the %LICENSE% license.

import logging

try:
    QString = unicode
except NameError:
    # Python 3
    QString = str
logger = logging.getLogger(__name__)
import argparse

# conda install -c conda-forge begins
# import begin

import PyQt5

from PyQt5.QtWidgets import (QGridLayout, QLabel, QPushButton, QLineEdit,
                         QApplication, QWidget, QGridLayout, QSpinBox,
                         QLineEdit, QCheckBox, QComboBox, QTextEdit, QDialog,
                         QMainWindow, QDoubleSpinBox)

from PyQt5 import QtGui, QtWidgets
import sys
import os.path as op
import copy
from .. import dictwidgetqt, iowidgetqt
from . import unconnected_cylinders
from imma import dili
def complicated_to_yaml(cfg):
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


class CylindersWidget(QtWidgets.QWidget):
    def __init__(self, ncols=2):
        super(CylindersWidget, self).__init__()
        self.ncols = ncols

        self.config = dili.get_default_args(unconnected_cylinders.UnconnectedCylinderGenerator)
        self.gen = None
        self.init_ui()

    def run(self):
        new_cfg = self.configwg.config_as_dict()
        logger.debug(str(new_cfg))
        self.config = new_cfg
        self.gen = cylinders.CylinderGenerator(**self.config)
        self.gen.run()

    def _show_stats(self):
        df = self.gen.get_stats()
        from .. import tablewidget


        dfmerne = df[["length", "volume", "surface"]].sum() / self.gen.area_volume
        dfmernef = dfmerne.to_frame().transpose()
        # dfmernef.insert(0, "", dfmernef.index)
        # import ipdb; ipdb.set_trace()
        tw = tablewidget.TableWidget(self, dataframe=dfmernef)
        self.mainLayout.addWidget(tw)

        from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

        # from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
        import matplotlib.pyplot as plt


        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        # self.toolbar = NavigationToolbar(self.canvas, self)
        self.mainLayout.addWidget(self.canvas)
        plt.subplot(141)
        df[["length"]].boxplot()
        plt.subplot(142)
        df[['radius']].boxplot()

        plt.subplot(143)
        df[["surface"]].boxplot()

        plt.subplot(144)
        df[["volume"]].boxplot()

        # TODO take care about redrawing
        dfdescribe = df.describe()
        dfdescribe.insert(0, "", dfdescribe.index)
        tw = tablewidget.TableWidget(self, dataframe=dfdescribe)
        tw.show()
        tw.raise_()
        tw.setMinimumWidth(600)
        tw.setMinimumHeight(200)

        self.mainLayout.addWidget(tw, 0, 2, 5, 2)
        self.resize(600, 700)

    def init_ui(self):
        self.mainLayout = QGridLayout(self)
        hide_keys = ["build", "gtree"]

        self.configwg = dictwidgetqt.DictWidget(self.config, hide_keys=hide_keys)
        self.mainLayout.addWidget(self.configwg)

        # self.mainLayout.setColumnMinimumWidth(text_col, 500)

        btn_accept = QPushButton("Run", self)
        btn_accept.clicked.connect(self.btnAccept)
        self.mainLayout.addWidget(btn_accept)  # , (gd_max_i / 2), text_col)

        self.ui_output_dir_widget = iowidgetqt.SetDirWidget("~", "output directory")
        self.mainLayout.addWidget(self.ui_output_dir_widget)  # , (gd_max_i / 2), text_col)

        btn_save = QPushButton("Save", self)
        btn_save.clicked.connect(self.btnSave)
        self.mainLayout.addWidget(btn_save)  # , (gd_max_i / 2), text_col)
        # self.config.updated.connect(self.on_config_update)

    def btnAccept(self):

        logger.debug("btnAccept")
        logger.debug(str(self.config))
        self.run()
        self._show_stats()

    def btnSave(self):
        # filename = "file{:05d}.jpg"
        filename = "file%05d.jpg"
        # filename = QtGui.QFileDialog.getSaveFileName(
        #     self,
        #     "Save file",
        #     init_filename,
        #     ""
        # )
        # filename = str(filename)

        if self.gen is None:
            self.run()
            self._show_stats()

        filename = op.join(self.ui_output_dir_widget.get_dir(), filename)
        filename = iowidgetqt.str_format_old_to_new(filename)
        self.gen.saveVolumeToFile(filename=filename)

    def on_config_update(self):
        pass

    def config_as_dict(self):
        dictionary = self.config.as_dict()
        for key, value in dictionary.items():
            from PyQt5.QtCore import pyqtRemoveInputHook

            pyqtRemoveInputHook()
            # import ipdb; ipdb.set_trace() #  noqa BREAKPOINT
            if type(value) == QString:
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
        # required=True,
        help='input file'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='Debug mode')
    args = parser.parse_args()

    if args.debug:
        ch.setLevel(logging.DEBUG)

    app = QApplication(sys.argv)
    cw = CylindersWidget()
    cw.show()
    app.exec_()


if __name__ == "__main__":
    main()
