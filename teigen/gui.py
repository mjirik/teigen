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
import inspect
import collections
import numpy as np
import scipy

import dictwidgetqt, iowidgetqt
import generators.cylinders
import generators.gensei_wrapper

from pyqtconfig import ConfigManager

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
        self.dataframes = {}
        self.figures = {}
        self.ui_stats_shown = False
        self.teigen = Teigen()
        self.version = "0.1.19"
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
        if self.ui_stats_shown:
        #     self._wg_tab_describe.deleteLater()
        #     self._wg_tab_describe = None
        #     self._wg_tab_merne.deleteLater()
        #     self._wg_tab_merne = None
            pass
        else:

            self.stats_tab_wg = QTabWidget()
            self.mainLayout.addWidget(self.stats_tab_wg, 0, 3, 5, 2)
        if True:
            from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
            # from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
            import matplotlib.pyplot as plt

            self.figure = plt.figure()
            self.canvas = FigureCanvas(self.figure)
            # self.toolbar = NavigationToolbar(self.canvas, self)
            self.stats_tab_wg.addTab(self.canvas, 'Graphs')



        df = self.gen.getStats()
        import tablewidget
        to_rename = {
            "length": "length [mm]",
            "volume": "volume [mm^3]",
            "surface": "surface [mm^2]",
            "radius": "radius [mm^2]"
        }
        to_rename_density = {
            "length": "length [mm^-2]",
            "volume": "volume []",
            "surface": "surface [mm^-1]",
            "radius": "radius [mm^-2]"
        }

        dfmerne = df[["length", "volume", "surface", "radius"]].sum() / self.gen.area_volume
        print "merne"
        print dfmerne
        dfmernef = dfmerne.to_frame().transpose().rename(columns=to_rename_density)
        # dfmernef.insert(0, "", dfmernef.index)
        # import ipdb; ipdb.set_trace()

        self._wg_tab_merne = tablewidget.TableWidget(self, dataframe=dfmernef)
        self.stats_tab_wg.addTab(self._wg_tab_merne, "Density table")
        self.dataframes["density"] = dfmernef

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
        self.dataframes["describe"] = dfdescribe

        self._wg_tab_describe = tablewidget.TableWidget(self, dataframe=dfdescribe)
        self._wg_tab_describe.show()
        self._wg_tab_describe.raise_()
        self._wg_tab_describe.setMinimumWidth(600)
        self._wg_tab_describe.setMinimumHeight(200)

        # self.mainLayout.addWidget(self._wg_tab_describe, 0, 2, 5, 2)
        self.stats_tab_wg.addTab(self._wg_tab_describe, "Stats table")
        self.resize(600,700)

        self.ui_stats_shown = True




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
        wtitle = "Teigen " + self.version
        self.setWindowTitle(wtitle)
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

        postprocessing_params = dictwidgetqt.get_default_args(self.teigen.postprocessing)
        self.posprocessing_wg = dictwidgetqt.DictWidget(postprocessing_params)
        self.mainLayout.addWidget(self.posprocessing_wg)


        btn_accept = QPushButton("Run", self)
        btn_accept.clicked.connect(self.btnRun)
        self.mainLayout.addWidget(btn_accept) # , (gd_max_i / 2), text_col)

        self.ui_output_dir_widget = iowidgetqt.SetDirWidget("~/teigen_data/slice{:06d}.jpg", "output directory")
        self.ui_output_dir_widget.setToolTip("Data are stored in defined directory.\nOutput format is based on file extension.\nFor saving into image stack use 'filename{:06d}.jpg'")
        self.mainLayout.addWidget(self.ui_output_dir_widget) # , (gd_max_i / 2), text_col)

        btn_save = QPushButton("Save", self)
        btn_save.setToolTip("Save image slices and meta information")
        btn_save.clicked.connect(self.btnSave)
        self.mainLayout.addWidget(btn_save) # , (gd_max_i / 2), text_col)
        # self.config.updated.connect(self.on_config_update)


        ## pyqtgraph experiments
        import dictwidgetpyqtgraph
        import pyqtgraph
        from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
        input_params = {
            "voxelsize": [0.01, 0.01, 0.01],
            'area_size_px': [100, 100, 100],
            'area_size': [10.0, 10.0, 10.0]
        }
        properties = {
            "children": {
                "voxelsize_mm": {
                    "title": 'voxelsize [mm]',
                    "children": {
                        "0": {
                            "title": "z",
                            'suffix': 'm',
                            'siPrefix': True
                        },
                        "1": {
                            "title": "x",
                            'suffix': 'm',
                            'siPrefix': True
                        },
                        "2": {
                            "title": "y",
                            'suffix': 'm',
                            'siPrefix': True
                        }
                    }
                }
            }
        }

        gr_struct = dictwidgetpyqtgraph.to_pyqtgraph_struct('params', input_params, properties=properties)

        gr_struct['children'].append(
            dictwidgetpyqtgraph.ComplexParameter(name='Custom parameter group (reciprocal values)'))
        p = Parameter.create(**gr_struct)

        # t = ParameterTree()
        # t.setParameters(p, showTop=False)
        # t.setMinimumWidth(300)
        # t.show()
        # self.mainLayout.addWidget(t, 0, 1, 5, 1)

        ## end of pyqtgraph tree
    def btnRun(self):

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

        # postprocessing
        if "generate_volume" in dir(self.gen):
            self.teigen.data3d = self.gen.generate_volume()
            self.teigen.voxelsize_mm = self.gen.voxelsize_mm
            postprocessing_params = self.posprocessing_wg.config_as_dict()
            data3d = self.teigen.postprocessing(**postprocessing_params)
            self.gen.data3d = data3d

        # filename = op.join(self.ui_output_dir_widget.get_dir(), filename)
        filename = self.ui_output_dir_widget.get_dir()

        filename = iowidgetqt.str_format_old_to_new(filename)



        self.gen.saveVolumeToFile(filename)
        fn_base = self.filename_base(filename)
        for dfname in self.dataframes:
            df = self.dataframes[dfname].to_csv(fn_base+"_" + dfname + ".csv")
        self.figure.savefig(fn_base + "_" + "graph.pdf")
        self.figure.savefig(fn_base + "_" + "graph.png")

    def filename_base(self, filename):

        import re
        filename = re.sub(r"({.*})", r"", filename)

        root, ext = op.splitext(filename)
        return root


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


class Teigen():
    def __init__(self):
        self.data3d = None
        self.voxelsize_mm = None


    def postprocessing(self, gaussian_filter=True, gaussian_filter_sigma_mm=1.0, gaussian_noise=True, gaussian_noise_stddev=10.0, gaussian_noise_center=0.0, limit_negative_intensities=True):
        if gaussian_filter:
            sigma_px = gaussian_filter_sigma_mm / self.voxelsize_mm

            self.data3d = scipy.ndimage.filters.gaussian_filter(
                self.data3d,
                sigma=sigma_px)

        if gaussian_noise:
            dt = self.data3d.dtype
            noise = np.random.normal(loc=gaussian_noise_center, scale=gaussian_noise_stddev, size=self.data3d.shape)
            self.data3d = (self.data3d + noise).astype(self.data3d.dtype)

        if limit_negative_intensities:
            self.data3d[self.data3d < 0] = 0

        return self.data3d

    # def save_volume_to_file(self, filename):
    #
    #     import io3
    #     import io3d.misc
    #     import numpy as np
    #     data = {
    #         'data3d': self.data3d.astype(np.uint8), #* self.output_intensity,
    #         'voxelsize_mm': self.voxelsize_mm,
    #         # 'segmentation': np.zeros_like(self.data3d, dtype=np.int8)
    #     }
    #     io3d.write(data, filename)

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