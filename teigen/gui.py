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
        QComboBox, QTextEdit, QDialog, QMainWindow, QDoubleSpinBox, QTabWidget, QStatusBar

from PyQt4 import QtGui
import sys
import os.path as op
import copy
import inspect
import collections
import numpy as np
import scipy
import re

import dictwidgetqt
import iowidgetqt
import dictwidgetpyqtgraph
import generators.cylinders
import generators.gensei_wrapper
import io3d.datawriter
import io3d.misc

from pyqtconfig import ConfigManager


class TeigenWidget(QtGui.QWidget):
    def __init__(self, ncols=2, qapp=None):
        super(TeigenWidget, self).__init__()
        self.ncols = ncols
        self.gen = None
        self.dataframes = {}
        self.figures = {}
        self.ui_stats_shown = False
        self.teigen = Teigen()
        self.version = self.teigen.version
        self.config = {}
        self.run_number = 0
        self.qapp = qapp
        self.init_ui()

    def collect_config(self):
        print "generator args"
        id = self.gen_tab_wg.currentIndex()
        config = self._ui_generator_widgets[id].config_as_dict()
        logger.debug(str(config))

        none, area_cfg = dictwidgetpyqtgraph.from_pyqtgraph_struct(self.area_sampling_params.saveState())

        config.update(area_cfg["Area Sampling"])
        filepattern = self.ui_output_dir_widget.get_dir()
        series_number = io3d.datawriter.get_unoccupied_series_number(filepattern=filepattern)
        config["filepattern"] = filepattern
        config['filepattern_series_number'] = series_number
        config["generator_id"] = id

        config["postprocessing"] = self.posprocessing_wg.config_as_dict()
        config["required_teigen_version"] = self.teigen.version
        self.config = config


    def run(self):

        self.collect_config()

        # self.config = new_cfg
        self.teigen.run(**self.config)



    def _show_stats(self):
        to_rename = {
            "length": "length [mm]",
            "volume": "volume [mm^3]",
            "surface": "surface [mm^2]",
            "radius": "radius [mm]"
        }
        to_rename_density = {
            "length": "length [mm^-2]",
            "volume": "volume []",
            "surface": "surface [mm^-1]"
            # "radius": "radius [mm^-2]"
        }

        run_number_alpha = chr(ord("A") + self.run_number)
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
            self.stats_tab_wg.addTab(self.canvas, 'Graphs ' + run_number_alpha)

        df = self.teigen.gen.getStats()

        plt.subplot(141)
        df[["length"]].rename(columns=to_rename).boxplot(return_type='axes')
        plt.subplot(142)
        df[['radius']].rename(columns=to_rename).boxplot(return_type='axes')

        plt.subplot(143)
        df[["surface"]].rename(columns=to_rename).boxplot(return_type='axes')

        plt.subplot(144)
        df[["volume"]].rename(columns=to_rename).boxplot(return_type='axes')
        self.figure.tight_layout()

        import tablewidget

        dfmerne = df[["length", "volume", "surface"]].sum() / self.teigen.gen.area_volume
        print "merne"
        print dfmerne
        dfmernef = dfmerne.to_frame().transpose().rename(columns=to_rename_density)
        # dfmernef.insert(0, "", dfmernef.index)
        # import ipdb; ipdb.set_trace()

        self._wg_tab_merne = tablewidget.TableWidget(self, dataframe=dfmernef)
        # self.stats_tab_wg.addTab(self._wg_tab_merne, "Density table")
        self.dataframes["density"] = dfmernef


        # TODO take care about redrawing
        dfdescribe = df.describe()
        dfdescribe.insert(0, "", dfdescribe.index)
        self.dataframes["describe"] = dfdescribe

        self._wg_tab_describe = tablewidget.TableWidget(self, dataframe=dfdescribe)
        self._wg_tab_describe.setMinimumWidth(800)
        self._wg_tab_describe.setMinimumHeight(200)

        self._wg_tables = QtGui.QWidget()
        self._wg_tables.setLayout(QGridLayout())
        self._wg_tables.layout().addWidget(self._wg_tab_describe)
        self._wg_tables.layout().addWidget(self._wg_tab_merne)

        self._wg_tab_describe.show()
        self._wg_tab_describe.raise_()
        # self.mainLayout.addWidget(self._wg_tab_describe, 0, 2, 5, 2)
        # self.stats_tab_wg.addTab(self._wg_tab_describe, "Stats table")
        self.stats_tab_wg.addTab(self._wg_tables, "Sumary " + run_number_alpha)
        # self.resize(600,700)



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

        self.statusBar = QStatusBar()
        self.mainLayout.addWidget(self.statusBar, 10, 0, 1, 2)
        self.progressBar = QtGui.QProgressBar()
        self.progressBar.setRange(0, 10000)
        self.progressBar.setValue(0)

        self.statusBar.addWidget(self.progressBar)
        self.progressBar.show()



        hide_keys = ["build", "gtree", "voxelsize_mm", "areasize_px", "resolution", "n_slice", "dims"]
        self.gen_tab_wg = QTabWidget()
        self.mainLayout.addWidget(self.gen_tab_wg, 0, 1)

        rename_captions_dict = {
            "voxelsize_mm": "voxel size [mm]",
            }

        self._ui_generator_widgets = []
        for i, config in enumerate(self.teigen.configs):
            wg = dictwidgetqt.DictWidget(
                self.teigen.configs[i],
                hide_keys=hide_keys,
                captions=rename_captions_dict,
                ncols=1,
            )
            self._ui_generator_widgets.append(wg)
            self.gen_tab_wg.addTab(wg, self.teigen.generators_names[i])
        # self.gen_tab_wg.addTab(gen_wg, "cylinder generator")
        # self.gen_tab_wg.addTab(gen_wg, "gensei generator")

        # self.mainLayout.setColumnMinimumWidth(text_col, 500)

        self.ui_output_dir_widget = iowidgetqt.SetDirWidget("~/teigen_data/{seriesn:03d}/slice{:06d}.jpg", "output directory")
        self.ui_output_dir_widget.setToolTip("Data are stored in defined directory.\nOutput format is based on file extension.\nFor saving into image stack use 'filename{:06d}.jpg'")
        self.mainLayout.addWidget(self.ui_output_dir_widget, 1, 1) # , (gd_max_i / 2), text_col)
        btn_accept = QPushButton("Generate skeleton", self)

        btn_accept.clicked.connect(self.btnRun)
        self.mainLayout.addWidget(btn_accept, 2, 1) # , (gd_max_i / 2), text_col)

        postprocessing_params = dictwidgetqt.get_default_args(self.teigen.postprocessing)
        self.posprocessing_wg = dictwidgetqt.DictWidget(postprocessing_params)
        self.mainLayout.addWidget(self.posprocessing_wg, 3, 1)




        btn_save = QPushButton("Generate volumetric data", self)
        btn_save.setToolTip("Save image slices and meta information")
        btn_save.clicked.connect(self.btnSave)
        self.mainLayout.addWidget(btn_save, 4, 1) # , (gd_max_i / 2), text_col)
        # self.config.updated.connect(self.on_config_update)


        ## pyqtgraph experiments
        import dictwidgetpyqtgraph
        import pyqtgraph
        from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
        input_params = {
            "Area Sampling":  dictwidgetpyqtgraph.AreaSamplingParameter(name='Area Sampling'),
            # TODO add more lines here
            # "Intensity Profile": dictwidgetpyqtgraph.ScalableFloatGroup(
            #     name="Intensity Profile", children=[
            #         {'name': '0.2', 'type': 'float', 'value': "100"},
            #         {'name': '0.4', 'type': 'float', 'value': "115"},
            #     ])
        }
        gr_struct = dictwidgetpyqtgraph.to_pyqtgraph_struct('params', input_params, opts={})
        p = Parameter.create(**gr_struct)

        t = ParameterTree()
        t.setParameters(p, showTop=False)
        t.setMinimumWidth(300)
        t.show()
        self.mainLayout.addWidget(t, 0, 0, 5, 1)
        self.area_sampling_wg = t
        self.area_sampling_params = p

        self.teigen.progress_callback = self._progressbar_update

    def _progressbar_update(self, obj, level, *args, **kwargs):
        self.progressBar.setValue(int(10000*level))
        ## end of pyqtgraph tree

    def btnRun(self):

        logger.debug("btnAccept")
        # logger.debug(str(self.config))
        self.run()
        self._show_stats()
        self.run_number += 1

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

        if self.teigen.need_run:
            self.run()
            self._show_stats()


        # filename = op.join(self.ui_output_dir_widget.get_dir(), filename)
        filename = self.config["filepattern"]

        # filename = iowidgetqt.str_format_old_to_new(filename)



        self.teigen.save_volume()
        fn_base, fn_ext = self.teigen.filepattern_split()
        for dfname in self.dataframes:
            df = self.dataframes[dfname].to_csv(fn_base+"_" + dfname + ".csv")
        self.figure.savefig(fn_base + "_" + "graph.pdf")
        self.figure.savefig(fn_base + "_" + "graph.png")
        self.figure.savefig(fn_base + "_" + "graph.svg")

        # self.teigen.gen.saveVolumeToFile(filename)



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
        self.version = "0.1.21"
        self.data3d = None
        self.voxelsize_mm = None
        self.need_run = True
        self.gen = None
        self.generators_classes =  [
            generators.cylinders.CylinderGenerator,
            generators.gensei_wrapper.GenseiGenerator
        ]
        self.generators_names = [
            "Cylinder generator",
            "Gensei generator"
        ]
        self.configs = [dictwidgetqt.get_default_args(conf) for conf in self.generators_classes]
        self.config = self.configs[0]
        self.progress_callback = None

    def run(self, **config):
        import io3d.misc
        self.config = copy.deepcopy(config)

        if "required_teigen_version" in config.keys():
            reqired_version = config["required_teigen_version"]
            if  reqired_version != self.version:
                logger.error(
                    "Wrong teigen version. Required: " + reqired_version + " , actual " + self.version)
                return

        filepattern = config["filepattern"]
        if "filepattern_series_number" in config.keys():
            series_number = config["filepattern_series_number"]
            config['filepattern'] = io3d.datawriter.filepattern_fill_series_number(
                filepattern,
                series_number=series_number
            )




        id = config.pop('generator_id')
        cfg_export_fcn = [
            self._area_sampling_general_export,
            self._area_sampling_gensei_export,
        ]

        area_dct = dictwidgetqt.subdict(config, ["voxelsize_mm", "areasize_mm", "areasize_px"])
        config.pop("voxelsize_mm")
        config.pop("areasize_mm")
        config.pop("areasize_px")

        cfg = cfg_export_fcn[id](area_dct)
        config.update(cfg)

        generator_class = self.generators_classes[id]
        # self.config = get_default_args(generator_class)

        # select only parameters for generator
        generator_default_config = dictwidgetqt.get_default_args(generator_class)
        generator_config = dictwidgetqt.subdict(config,generator_default_config.keys())
        self.gen = generator_class(**generator_config)
        self.gen.progress_callback = self.progress_callback
        # self.gen = generators.gensei_wrapper.GenseiGenerator(**self.config2)
        # self.gen = generators.gensei_wrapper.GenseiGenerator()
        self.gen.run()

        self.__generate_vtk()

        self.need_run = False

    def __generate_vtk(self, vtk_file="~/tree.vtk"):
        vtk_file = op.expanduser(vtk_file)
        from tree import TreeBuilder

        if "tree_data" in dir(self.gen):

            tvg = TreeBuilder('vtk')
            # yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
            # tvg.importFromYaml(yaml_path)
            tvg.voxelsize_mm = self.voxelsize_mm
            tvg.shape = self.gen.areasize_px
            tvg.tree_data = self.gen.tree_data
            output = tvg.buildTree() # noqa
            # tvg.show()
            tvg.saveToFile(vtk_file)

    def filepattern_split(self):
        """
        Return base and ext of file. The slice_number and slice_position is ignored.
        :return:
        """
        import io3d.datawriter
        filepattern = self.config["filepattern"]
        filepattern_series_number = self.config["filepattern_series_number"]

        # filepattern = re.sub(r"({\s*slicen\s*:?.*})", r"", filepattern)
        # filepattern = re.sub(r"({\s*slice_number\s*:?.*})", r"", filepattern)
        # filepattern = re.sub(r"({\s*slicep\s*:?.*})", r"", filepattern)
        # filepattern = re.sub(r"({\s*slice_position\s*:?.*})", r"", filepattern)

        filepattern = re.sub(r"({\s*})", r"", filepattern)

        filepattern = io3d.datawriter.filepattern_fill_series_number(filepattern, filepattern_series_number)
        filepattern = re.sub(r"({.*?})", r"", filepattern)
        root, ext = op.splitext(filepattern)
        return root, ext

    def save_volume(self):
        import io3d.misc
        fn_base, fn_ext = self.filepattern_split()
        io3d.misc.obj_to_file(self.config, filename=fn_base + "_parameters.yaml")
        # postprocessing
        if "generate_volume" in dir(self.gen):
            self.data3d = self.gen.generate_volume()
            self.voxelsize_mm = self.gen.voxelsize_mm
            postprocessing_params = self.config["postprocessing"]
            data3d = self.postprocessing(**postprocessing_params)
            self.gen.data3d = data3d
        self.gen.saveVolumeToFile(self.config["filepattern"])

    def postprocessing(
            self,
            gaussian_blur=True,
            gaussian_filter_sigma_mm=1.0,
            add_noise=True,
            # gaussian_noise_stddev=10.0,
            # gaussian_noise_center=0.0,
            limit_negative_intensities=True,
            noise_random_generator_seed=0,
            exponent=-1,
            lambda_start=0,
            lambda_range=-1,
            noise_amplitude = 40.0,
            noise_mean = 30.0
    ):
        if gaussian_blur:
            sigma_px = gaussian_filter_sigma_mm / self.voxelsize_mm

            self.data3d = scipy.ndimage.filters.gaussian_filter(
                self.data3d,
                sigma=sigma_px)

        if add_noise:
            import ndnoise.generator

            dt = self.data3d.dtype
            noise = ndnoise.noises(
                shape=self.data3d.shape,
                sample_spacing=self.voxelsize_mm,
                exponent=exponent,
                random_generator_seed=noise_random_generator_seed,
                lambda_start=lambda_start,
                lambda_range=lambda_range

            ).astype(np.float16)
            mx = np.max(noise)
            noise = noise_amplitude * noise/mx
            noise += noise_mean
            noise = noise.astype(self.data3d.dtype)
            # noise = np.random.normal(loc=gaussian_noise_center, scale=gaussian_noise_stddev, size=self.data3d.shape)
            self.data3d = (self.data3d + noise).astype(self.data3d.dtype)

        if limit_negative_intensities:
            self.data3d[self.data3d < 0] = 0

        return self.data3d

    def _area_sampling_gensei_export(self, area_sampling_params):
        asp = area_sampling_params
        vs_mm = np.asarray(asp["voxelsize_mm"])
        resolution = 1.0 / vs_mm
        dct = {
            'dims': asp["areasize_mm"],
            'n_slice': asp["areasize_px"][0],
            'resolution': [resolution[1], resolution[2]]
        }
        return dct


    def _area_sampling_general_export(self, area_sampling_params):
        return {
            'voxelsize_mm': area_sampling_params["voxelsize_mm"],
            'areasize_px': area_sampling_params["areasize_px"]
        }
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
        '-p', '--parameterfile',
        default=None,
        # required=True,
        help='input parameter file'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='Debug mode')
    args = parser.parse_args()

    if args.debug:
        ch.setLevel(logging.DEBUG)


    if args.parameterfile is None:
        app = QApplication(sys.argv)
        cw = TeigenWidget()
        cw.show()
        app.exec_()
    else:
        params = io3d.misc.obj_from_file(args.parameterfile)
        tg = Teigen()
        tg.run(**params)
        tg.save_volume()



if __name__ == "__main__":
    main()