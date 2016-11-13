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
# conda install -c conda-forge begins
# import begin

import PyQt4
from PyQt4.QtGui import QGridLayout, QLabel,\
    QPushButton, QLineEdit, QApplication, QWidget, QGridLayout, QSpinBox, QLineEdit, QCheckBox,\
        QComboBox, QTextEdit, QDialog, QMainWindow, QDoubleSpinBox

from PyQt4 import QtGui
import sys
import os.path
import copy

from .. import dictwidgetqt
import cylinders

from pyqtconfig import ConfigManager
import inspect
import collections

def get_default_args(obj):
    argspec = inspect.getargspec(obj.__init__)
    args = argspec.args[1:]
    defaults = argspec.defaults
    print "---- args"
    print args
    print defaults
    # import ipdb; ipdb.set_trace() #  noqa BREAKPOINT
    # dc = dict(zip(args, defaults))
    dc = collections.OrderedDict(zip(args, defaults))
    return dc



class CylindersWidget(QtGui.QWidget):
    def __init__(self, config_in, ncols=2):
        super(CylindersWidget, self).__init__()
        self.config_in = config_in
        self.ncols = ncols

        self.config = get_default_args(cylinders.CylinderGenerator)
        print "default args"
        print self.config
        self.init_ui()


    def run(self):
        print "generator args"
        print self.config
        self.gen = cylinders.CylinderGenerator(self.config)
        self.gen.run()

    def _show_stats(self):
        df = self.gen.getStats()
        dfdescribe = df.describe()
        # TODO take care about redrawing
        from .. import tablewidget
        tw = tablewidget.TableWidget(self, dataframe=dfdescribe)

        self.mainLayout.addWidget(tw)
        tw.show()
        tw.raise_()

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

        self.resize(600,700)




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

        self.configwg = dictwidgetqt.DictWidget(self.config)
        self.mainLayout.addWidget(self.configwg)

        # self.mainLayout.setColumnMinimumWidth(text_col, 500)

        btn_accept = QPushButton("Run", self)
        btn_accept.clicked.connect(self.btnAccept)
        self.mainLayout.addWidget(btn_accept) # , (gd_max_i / 2), text_col)

        # self.config.updated.connect(self.on_config_update)

    def btnAccept(self):

        logger.debug("btnAccept")
        logger.debug(str(self.config))
        logger.debug(str(self.configwg.config_as_dict()))
        self.run()
        self._show_stats()
        filename = QtGui.QFileDialog.getSaveFileName(
            self,
            "Save file",
            "file{:05d}.jpg",
            # ".jpg"
        )
        filename = str(filename)
        self.gen.saveVolumeToFile(filename=filename)



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
    cw = CylindersWidget(cfg)
    cw.show()
    app.exec_()


if __name__ == "__main__":
    main()