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
        QComboBox, QTextEdit, QDialog, QMainWindow, QDoubleSpinBox, QTabWidget

from PyQt4 import QtGui
import sys
import os.path as op
import copy

import dictwidgetqt, iowidgetqt
import generators.cylinders
import generators.gensei_wrapper

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


class TeigenWidget(QtGui.QWidget):
    def __init__(self, ncols=2):
        super(TeigenWidget, self).__init__()
        self.ncols = ncols

        self.generators_classes =  [
            generators.cylinders.CylinderGenerator,
            generators.gensei_wrapper.GenseiGenerator
        ]
        self.generators_names = [
            "Cylinder generator",
            "Gensei generator"
        ]
        self.configs = [get_default_args(conf) for conf in self.generators_classes]
        self.config = self.configs[0]

        print "default args"
        print self.config
        self.gen = None
        self.init_ui()


    def run(self):
        print "generator args"
        id = self.gen_tab_wg.currentIndex()
        new_cfg = self._ui_generator_widgets[id].config_as_dict()
        logger.debug(str(new_cfg))
        self.config = new_cfg
        generator_class = self.generators_classes[id]
        # self.config = get_default_args(generator_class)
        self.gen = generator_class(**self.config)
        # self.gen = generators.gensei_wrapper.GenseiGenerator(**self.config2)
        # self.gen = generators.gensei_wrapper.GenseiGenerator()
        self.gen.run()

    def _show_stats(self):
        df = self.gen.getStats()
        import tablewidget
        to_rename = {
            "length": "length [mm]",
            "volume": "volume [mm^3]",
            "surface": "surface [mm^2]",
            "radius": "radius [mm^2]"
        }
        to_rename_relative = {
            "length": "length [mm^-2]",
            "volume": "volume []",
            "surface": "surface [mm^-1]",
            "radius": "radius [mm^-2]"
        }

        dfmerne = df[["length", "volume", "surface", "radius"]].sum() / self.gen.area_volume
        print "merne"
        print dfmerne
        dfmernef = dfmerne.to_frame().transpose().rename(columns=to_rename_relative)
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
        df[["length"]].rename(columns=to_rename).boxplot(return_type='axes')
        plt.subplot(142)
        df[['radius']].rename(columns=to_rename).boxplot(return_type='axes')

        plt.subplot(143)
        df[["surface"]].rename(columns=to_rename).boxplot(return_type='axes')

        plt.subplot(144)
        df[["volume"]].rename(columns=to_rename).boxplot(return_type='axes')

        # TODO take care about redrawing
        dfdescribe = df.describe()
        dfdescribe.insert(0, "", dfdescribe.index)
        tw = tablewidget.TableWidget(self, dataframe=dfdescribe)
        tw.show()
        tw.raise_()
        tw.setMinimumWidth(600)
        tw.setMinimumHeight(200)

        self.mainLayout.addWidget(tw,0,2,5,2)
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

    def _get_generator(self, id):
        pass


    def init_ui(self):
        self.mainLayout = QGridLayout(self)
        hide_keys = ["build", "gtree"]
        self.gen_tab_wg = QTabWidget()
        self.mainLayout.addWidget(self.gen_tab_wg)

        rename_captions_dict = {
            "voxelsize_mm": "voxel size [mm]",
            }

        self._ui_generator_widgets = []
        for i, config in enumerate(self.configs):
            wg = dictwidgetqt.DictWidget(
                self.configs[i],
                hide_keys=hide_keys,
                captions=rename_captions_dict)
            self._ui_generator_widgets.append(wg)
            self.gen_tab_wg.addTab(wg, self.generators_names[i])
        # self.gen_tab_wg.addTab(gen_wg, "cylinder generator")
        # self.gen_tab_wg.addTab(gen_wg, "gensei generator")

        # self.mainLayout.setColumnMinimumWidth(text_col, 500)

        btn_accept = QPushButton("Run", self)
        btn_accept.clicked.connect(self.btnAccept)
        self.mainLayout.addWidget(btn_accept) # , (gd_max_i / 2), text_col)

        self.ui_output_dir_widget = iowidgetqt.SetDirWidget("~/teigen_data/slice{:06d}.jpg", "output directory")
        self.ui_output_dir_widget.setToolTip("Data are stored in defined directory.\nOutput format is based on file extension.\nFor saving into image stack use 'filename{:06d}.jpg'")
        self.mainLayout.addWidget(self.ui_output_dir_widget) # , (gd_max_i / 2), text_col)

        btn_save = QPushButton("Save", self)
        btn_save.setToolTip("Save image slices and meta information")
        btn_save.clicked.connect(self.btnSave)
        self.mainLayout.addWidget(btn_save) # , (gd_max_i / 2), text_col)
        # self.config.updated.connect(self.on_config_update)

    def btnAccept(self):

        logger.debug("btnAccept")
        logger.debug(str(self.config))
        self.run()
        self._show_stats()

    def btnSave(self):
        # filename = "file{:05d}.jpg"
        # filename = "file%05d.jpg"
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

        # filename = op.join(self.ui_output_dir_widget.get_dir(), filename)
        filename = self.ui_output_dir_widget.get_dir()

        filename = iowidgetqt.str_format_old_to_new(filename)



        self.gen.saveVolumeToFile(filename)



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
    cw = TeigenWidget()
    cw.show()
    app.exec_()


if __name__ == "__main__":
    main()