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

import logging.handlers
import PyQt4
from PyQt4.QtGui import QGridLayout, QLabel,\
    QPushButton, QLineEdit, QApplication, QWidget, QGridLayout, QSpinBox, QLineEdit, QCheckBox,\
        QComboBox, QTextEdit, QDialog, QMainWindow, QDoubleSpinBox, QTabWidget, QStatusBar

from PyQt4 import QtGui
import sys
import os
import os.path as op
import copy
import inspect
import collections
import numpy as np
import scipy
import re

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
# from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import matplotlib.pyplot as plt

import dictwidgetqt
import iowidgetqt
import dictwidgetpg
import generators.cylinders
import generators.gensei_wrapper
import generators.unconnected_cylinders
import io3d.datawriter
import io3d.misc
import dili

from pyqtconfig import ConfigManager


class TeigenWidget(QtGui.QWidget):
    def __init__(self, ncols=2, qapp=None, logfile="~/teigen.log", config=None):
        super(TeigenWidget, self).__init__()
        self.logfile = logfile
        self.ncols = ncols
        self.gen = None
        self.figures = {}
        self.ui_stats_shown = False
        self.teigen = Teigen(logfile=self.logfile)
        if config is not None:
            self.teigen.update_config(**config)
        self.version = self.teigen.version
        self.config = {}
        self.run_number = 0
        self.qapp = qapp
        self.init_ui()

    def collect_config_from_gui(self):
        print "generator args"
        id = self.gen_tab_wg.currentIndex()
        config = collections.OrderedDict()
        config["generators"] = collections.OrderedDict()
        for i, wg in enumerate(self._ui_generator_widgets):
            config["generators"][self.teigen.generators_names[i]] = wg.config_as_dict()

        # config = self._ui_generator_widgets[id].config_as_dict()
        logger.debug(str(config))

        none, area_cfg = dictwidgetpg.from_pyqtgraph_struct(self.area_sampling_params.saveState())

        # config.update(area_cfg["Area Sampling"])
        config["areasampling"] = area_cfg["Area Sampling"]
        filepattern = self.ui_output_dir_widget.get_dir()
        # series_number = io3d.datawriter.get_unoccupied_series_number(filepattern=filepattern)
        config["filepattern"] = filepattern
        # config['filepattern_series_number'] = series_number
        config["generator_id"] = id

        # config["postprocessing"] = self.posprocessing_wg.config_as_dict()
        config["postprocessing"] = area_cfg["Postprocessing"]
        config["required_teigen_version"] = self.teigen.version
        config["appearance"] = area_cfg["Appearance"]
        self.config = config

    def _parameters_changed(self):
        self.teigen.parameters_changed_before_save = True

    def collect_config_from_gui_and_push_to_teigen(self):
        self.collect_config_from_gui()
        self.teigen.update_config(**self.config)

    def run(self):
        self.collect_config_from_gui_and_push_to_teigen()


        # self.config = new_cfg
        self.teigen.run()

    def _show_stats(self):
        to_rename = {
            "length": "length [mm]",
            "volume": "volume [mm^3]",
            "surface": "surface [mm^2]",
            "radius": "radius [mm]"
        }
        to_rename_density = {
            "length": "length d. [mm^-2]",
            "volume": "volume d. []",
            "surface": "surface d. [mm^-1]"
            # "radius": "radius [mm^-2]"
        }

        # to_rename_density = {
        #     "length": "length [mm]",
        #     "volume": "volume [mm^3]",
        #     "surface": "surface [mm^2]"
        #     # "radius": "radius [mm^-2]"
        # }

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

        self.actual_subtab_wg = QTabWidget()
        self.stats_tab_wg.addTab(self.actual_subtab_wg, '' + run_number_alpha)
        if True:

            self.figure = plt.figure()
            self.canvas = FigureCanvas(self.figure)
            # self.toolbar = NavigationToolbar(self.canvas, self)
            self.actual_subtab_wg.addTab(self.canvas, 'Graphs ' + run_number_alpha)

        # df = self.teigen.gen.getStats()
        df = self.teigen.dataframes["elements"]

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

        # dfmernef.insert(0, "", dfmernef.index)
        # import ipdb; ipdb.set_trace()

        # TODO take care about redrawing
        # self.stats_tab_wg.addTab(self._wg_tab_merne, "Density table")

        dfdescribe = self.teigen.dataframes["describe"]
        dfmerne = self.teigen.dataframes["density"]
        dfoverall = self.teigen.dataframes["overall"]

        self._wg_tab_describe = tablewidget.TableWidget(self, dataframe=dfdescribe)
        self._wg_tab_merne = tablewidget.TableWidget(self, dataframe=dfmerne)
        self._wg_tab_overall = tablewidget.TableWidget(self, dataframe=dfoverall)
        self._wg_tab_describe.setMinimumWidth(800)
        self._wg_tab_describe.setMinimumHeight(200)
        self._wg_tab_merne.setMaximumHeight(80)
        self._wg_tab_overall.setMaximumHeight(80)

        self._wg_tables = QtGui.QWidget()
        self._wg_tables.setLayout(QGridLayout())
        self._wg_tables.layout().addWidget(self._wg_tab_describe)
        self._wg_tables.layout().addWidget(self._wg_tab_merne)
        self._wg_tables.layout().addWidget(self._wg_tab_overall)

        self._wg_tab_describe.show()
        self._wg_tab_describe.raise_()
        # self.mainLayout.addWidget(self._wg_tab_describe, 0, 2, 5, 2)
        # self.stats_tab_wg.addTab(self._wg_tab_describe, "Stats table")
        self.actual_subtab_wg.addTab(self._wg_tables, "Summary " + run_number_alpha)
        # self.resize(600,700)

        if self.teigen.polydata is not None:
            import imtools.show_segmentation_qt
            self._wg_show_3d = imtools.show_segmentation_qt.ShowSegmentationWidget(None, show_load_button=False)

            # self._wg_show_3d.add_vtk_file(op.expanduser(self.teigen.temp_vtk_file))
            self._wg_show_3d.add_vtk_polydata(self.teigen.polydata)
            self.actual_subtab_wg.addTab(self._wg_show_3d, "Visualization " + run_number_alpha)

        self.ui_stats_shown = True

        if (
                    self.teigen.config["postprocessing"]["noise_preview"] and
                    self.teigen.config["postprocessing"]["add_noise"]):
            self._noise_figure = plt.figure()
            self._noise_canvas = FigureCanvas(self._noise_figure)
            # self.toolbar = NavigationToolbar(self.canvas, self)
            self.actual_subtab_wg.addTab(self._noise_canvas, 'Noise ' + run_number_alpha)
            noise = self.teigen.generate_noise()
            plt.imshow(noise[0, :, :], cmap="gray")
            plt.colorbar()

    def update_stats(self):
        """
        Function is used after volume generation and save.
        :return:
        """
        import tablewidget
        self._wg_tab_overall.deleteLater()
        self._wg_tab_overall = None
        #     self._wg_tab_describe.deleteLater()
        # self._wg_tables
        dfoverall = self.teigen.dataframes["overall"]
        self._wg_tab_overall = tablewidget.TableWidget(self, dataframe=dfoverall)
        self._wg_tab_overall.setMaximumHeight(80)
        self._wg_tables.layout().addWidget(self._wg_tab_overall)
        #     self._wg_tab_describe.deleteLater()
        #     self._wg_tab_describe = None
        #     self._wg_tab_merne.deleteLater()
        #     self._wg_tab_merne = None

        # Show surface
        measurement_multiplier = self.teigen.config["postprocessing"]["measurement_multiplier"]
        surface_measurement = self.teigen.config["postprocessing"]["surface_measurement"]
        show_surface = self.teigen.config["appearance"]["show_surface"]
        if surface_measurement and measurement_multiplier > 0 and show_surface:
            fig = plt.figure()
            self._surface_figure = fig
            self._surface_canvas = FigureCanvas(self._surface_figure)
            # self.toolbar = NavigationToolbar(self.canvas, self)
            self.actual_subtab_wg.addTab(self._surface_canvas, 'Exact surface')

            # import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
            vertices = self.teigen._surface_vertices
            faces = self.teigen._surface_faces
            # fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111, projection='3d')

            # Fancy indexing: `verts[faces]` to generate a collection of triangles
            mesh = Poly3DCollection(vertices[faces])
            mesh.set_edgecolor('k')
            ax.add_collection3d(mesh)
            sh = self.teigen._numeric_surface_measurement_shape
            ax.set_xlim(0, sh[0])  # a = 6 (times two for 2nd ellipsoid)
            ax.set_ylim(0, sh[1])  # b = 10
            ax.set_zlim(0, sh[2])  # c = 16
            # plt.show()

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

        self.ui_stop_button = QPushButton("Stop", self)
        self.ui_stop_button.clicked.connect(self.btnStop)

        self.statusBar.addWidget(self.progressBar)
        self.statusBar.addWidget(self.ui_stop_button)
        self.progressBar.show()



        hide_keys = ["build", "gtree", "voxelsize_mm", "areasize_px", "resolution", "n_slice", "dims"]
        self.gen_tab_wg = QTabWidget()
        self.mainLayout.addWidget(self.gen_tab_wg, 0, 1, 1, 2)

        rename_captions_dict = {
            "voxelsize_mm": "voxel size [mm]",
            }

        # list is pointer. It causes problems with temporary reconstruction information
        # copy fix this issue
        teigen_config = copy.deepcopy(self.teigen.config)

        self._ui_generator_widgets = []
        for generator_name in teigen_config["generators"]:
            wg = dictwidgetqt.DictWidget(
                teigen_config["generators"][generator_name],
                hide_keys=hide_keys,
                captions=rename_captions_dict,
                ncols=1,
            )
            self._ui_generator_widgets.append(wg)
            self.gen_tab_wg.addTab(wg, generator_name)
        self.gen_tab_wg.setCurrentIndex(teigen_config["generator_id"])
        # self.mainLayout.setColumnMinimumWidth(text_col, 500)


        self.ui_output_dir_widget = iowidgetqt.SetDirWidget(
            self.teigen.config["filepattern"], "output directory")
        self.ui_output_dir_widget.setToolTip("Data are stored in defined directory.\nOutput format is based on file extension.\nFor saving into image stack use 'filename{:06d}.jpg'")
        self.mainLayout.addWidget(self.ui_output_dir_widget, 1, 1, 1, 2) # , (gd_max_i / 2), text_col)
        btn_accept = QPushButton("Step 1 - Preview - Generate skeleton", self)

        btn_save = QPushButton("Save parameters", self)
        btn_save.setToolTip("Save generator parameters")
        btn_save.clicked.connect(self.save_parameters)
        self.mainLayout.addWidget(btn_save, 2, 1, 1, 1) # , (gd_max_i / 2), text_col)

        btn_save = QPushButton("Save parameters and add to batch", self)
        btn_save.setToolTip("Save generator parameters and then add to batch")
        btn_save.clicked.connect(self.save_parameters_and_add_to_batch)
        self.mainLayout.addWidget(btn_save, 2, 2, 1, 1) # , (gd_max_i / 2), text_col)
        btn_accept.clicked.connect(self.btnRun)

        self.mainLayout.addWidget(btn_accept, 3, 1, 1, 2) # , (gd_max_i / 2), text_col)

        postprocessing_params = self.teigen.config["postprocessing"]
        # self.posprocessing_wg = dictwidgetqt.DictWidget(postprocessing_params)
        # self.mainLayout.addWidget(self.posprocessing_wg, 3, 1)

        btn_save = QPushButton("Step 2 - Generate and save volumetric data", self)
        btn_save.setToolTip("Save image slices and meta information")
        btn_save.clicked.connect(self.btnSave)
        self.mainLayout.addWidget(btn_save, 4, 1, 1, 2) # , (gd_max_i / 2), text_col)


        import pyqtgraph as pg
        ## pyqtgraph experiments
        import dictwidgetpg
        import pyqtgraph
        from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
        input_params = {
            "Area Sampling":  dictwidgetpg.AreaSamplingParameter(name='Area Sampling', **self.teigen.config["areasampling"]),
            "Postprocessing": postprocessing_params,
            "Batch processing": dictwidgetpg.BatchFileProcessingParameter(
                name="Batch processing" , children=[
                    {'name': 'Run batch', 'type': 'action'},
                ]),
            "Appearance": self.teigen.config["appearance"]
            # 'name': {'type': 'action'},
            # "dur": i5,
            # TODO add more lines here
            # "Intensity Profile": dictwidgetpyqtgraph.ScalableFloatGroup(
            #     name="Intensity Profile", children=[
            #         {'name': '0.2', 'type': 'float', 'value': "100"},
            #         {'name': '0.4', 'type': 'float', 'value': "115"},
            #     ])
        }
        gr_struct = dictwidgetpg.to_pyqtgraph_struct('params', input_params, opts={})
        # gr_struct["children"].append(i5)
        p = Parameter.create(**gr_struct)

        t = ParameterTree()
        t.setParameters(p, showTop=False)
        t.setMinimumWidth(350)
        # t.setColumnCount(3)
        t.show()

        p.sigStateChanged.connect(self._parameters_changed)
        p.param('Batch processing', 'Run batch').sigActivated.connect(self.run_batch)

        # how to add button
        # i5  = pg.TreeWidgetItem(["Item 5"])
        # b5 = QtGui.QPushButton('Button')
        # i5.setWidget(1, b5)
        # t.addTopLevelItem(i5)

        self.mainLayout.addWidget(t, 0, 0, 5, 1)
        self.area_sampling_wg = t
        self.area_sampling_params = p

        self.teigen.progress_callback = self._progressbar_update

    def _progressbar_update(self, obj, level, *args, **kwargs):
        self.progressBar.setValue(int(10000*level))
        if "statusbar_text" in kwargs:
            # add this in gui
            print "statusbar_text " + kwargs["statusbar_text"]
        ## end of pyqtgraph tree

    def run_batch(self):

        none, area_cfg = dictwidgetpg.from_pyqtgraph_struct(self.area_sampling_params.saveState())
        lst = area_cfg["Batch processing"].values()
        self.teigen.run_batch(lst)

    def save_parameters(self):
        self.collect_config_from_gui_and_push_to_teigen()
        self.teigen.save_parameters()

    def save_parameters_and_add_to_batch(self):
        self.collect_config_from_gui_and_push_to_teigen()
        self.teigen.save_parameters()
        fn_base = self.teigen._get_fn_base()

        config_filename = fn_base + "_parameters.yaml"
        self.area_sampling_params.param("Batch processing").add_filename(config_filename)

    def btnRun(self):

        logger.debug("btnAccept")
        # logger.debug(str(self.config))
        self.run()
        self._show_stats()
        self.run_number += 1

    def btnStop(self):
        pass

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
        self.figure.savefig(fn_base + "_" + "graph.pdf")
        self.figure.savefig(fn_base + "_" + "graph.png")
        self.figure.savefig(fn_base + "_" + "graph.svg")

        # self.teigen.gen.saveVolumeToFile(filename)
        self.update_stats()



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
    def __init__(self, logfile='~/tegen.log', loglevel=logging.DEBUG):
        self.config_file_manager = ConfigFileManager("teigen")
        self.config_file_manager.init_config_dir()
        self.loglevel = loglevel

        logger = logging.getLogger()
        handler = logging.handlers.RotatingFileHandler(
            op.expanduser(logfile),
            maxBytes=1000000,
            backupCount=9
        )
        handler.setLevel(self.loglevel)
        # formatter = logging.Formatter('%(asctime)s %(name)-18s %(levelname)-8s %(message)s')
        self.formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(name)-18s %(lineno)-5d %(funcName)-12s %(message)s')
        handler.setFormatter(self.formatter)
        logger.addHandler(handler)

        # streamhandler = logging.StreamHandler()
        # streamhandler.setFormatter(formatter)
        # self.memoryhandler = logging.handlers.MemoryHandler(1024*10, logging.DEBUG, streamhandler)
        self.memoryhandler = logging.handlers.MemoryHandler(1024*100)#, logging.DEBUG, streamhandler)
        self.memoryhandler.setLevel(self.loglevel)
        logger.addHandler(self.memoryhandler)

        logger.info("Starting Teigen")



        self.logfile=logfile
        self.version = "0.2.7"
        self.data3d = None
        self.voxelsize_mm = None
        self.need_run = True
        self.gen = None
        self.generators_classes =  [
            generators.cylinders.CylinderGenerator,
            generators.gensei_wrapper.GenseiGenerator,
            generators.cylinders.CylinderGenerator,
            generators.unconnected_cylinders.UnconnectedCylinderGenerator,
        ]
        self.generators_names = [
            "Cylinder generator",
            "Gensei generator",
            "Cylinder continues",
            "Unconnected cylinders"
        ]
        self._cfg_export_fcn = [
            self._area_sampling_general_export,
            self._area_sampling_gensei_export,
            self._area_sampling_general_export,
            self._area_sampling_general_export,
        ]
        self.use_default_config()
        self.progress_callback = None
        self.temp_vtk_file = op.expanduser("~/tree.vtk")
        # 3D visualization data, works for some generators
        self.polydata = None
        self.dataframes = {}
        self.stats_times = {}
        self.parameters_changed_before_save = True

    def use_default_config(self):
        self.config = self.get_default_config()

    def get_default_config(self):

        config = {}
        # self.config["generators"] = [dictwidgetqt.get_default_args(conf) for conf in self.generators_classes]

        hide_keys = ["build", "gtree", "voxelsize_mm", "areasize_px", "resolution",
                     "n_slice", "dims"]
        config["generators"] = collections.OrderedDict()
        for generator_cl, generator_name in zip(
                self.generators_classes,
                self.generators_names
        ):
            generator_params = dictwidgetqt.get_default_args(generator_cl)
            generator_params = dili.kick_from_dict(generator_params, hide_keys)
            config["generators"][generator_name] = generator_params

        # self.config["generator_id"] = self.generators_names[0]
        config["generator_id"] = 0
        # self.config = self.configs[0]
        config["postprocessing"] = dictwidgetqt.get_default_args(self.postprocessing)
        # config["postprocessing"][""] = dictwidgetqt.get_default_args(self.postprocessing)
        config["areasampling"] = {
            "voxelsize_mm": [1.0, 1.0, 1.0],
            "areasize_mm": [110.0, 100.0, 100.0],
            "areasize_px": [110, 100, 100]
        }
        config["filepattern"] = "~/teigen_data/{seriesn:03d}/slice{:06d}.jpg"
        # config['filepattern_series_number'] = series_number
        # self.config["voxelsize_mm"] = [1.0, 1.0, 1.0]
        # self.config["areasize_mm"] = [100.0, 100.0, 100.0]
        # self.config["areasize_px"] = [100, 100, 100]
        config["appearance"] = {
            "show_surface": True
        }
        return config

    def update_config(self, **config):
        import io3d.misc

        if "required_teigen_version" in config.keys():
            reqired_version = config["required_teigen_version"]
            if  reqired_version != self.version:
                logger.error(
                    "Wrong teigen version. Required: " + reqired_version + " , actual " + self.version)
                return
        config = copy.deepcopy(config)
        # there can be stored more than our config. F.e. some GUI dict reconstruction information
        self.config = dili.recursive_update(self.config, config)
        self.parameters_changed_before_save = True

    def run(self):
        import time

        t0 = time.time()
        config = copy.deepcopy(self.config)
        # filepattern = config["filepattern"]
        # if "filepattern_series_number" in config.keys():
        #     series_number = config["filepattern_series_number"]
        #     config['filepattern'] = io3d.datawriter.filepattern_fill_series_number(
        #         filepattern,
        #         series_number=series_number
        #     )
        self.config_file_manager.save_init(self.config)




        id = config.pop('generator_id')

        # area_dct = dili.subdict(config, ["voxelsize_mm", "areasize_mm", "areasize_px"])
        # config.pop("voxelsize_mm")
        # config.pop("areasize_mm")
        # config.pop("areasize_px")

        area_dct = config["areasampling"]

        area_cfg = self._cfg_export_fcn[id](area_dct)

        # TODO probably unused
        config.update(area_cfg)

        generator_class = self.generators_classes[id]
        # self.config = get_default_args(generator_class)

        # select only parameters for generator
        # generator_default_config = dictwidgetqt.get_default_args(generator_class)
        # generator_config = dictwidgetqt.subdict(config["generators"][id], generator_default_config.keys())
        generator_config = config["generators"].items()[id][1]
        generator_config.update(area_cfg)
        self.gen = generator_class(**generator_config)
        if id == 2:
            self.gen.MAKE_IT_SHORTER_CONSTANT = 0.0
            self.gen.OVERLAPS_ALOWED = True
        self.gen.progress_callback = self.progress_callback
        # self.gen = generators.gensei_wrapper.GenseiGenerator(**self.config2)
        # self.gen = generators.gensei_wrapper.GenseiGenerator()
        self.gen.run()

        self.polydata = self.__generate_vtk(self.temp_vtk_file)

        self.prepare_stats()


        t1 = time.time()
        self.time_run = t1 - t0
        logger.info("time: " + str(self.time_run))
        self.need_run = False
        self.parameters_changed_before_save = False

    def get_config_file_pattern(self):
        filepattern = self.config["filepattern"]
        filepattern = io3d.datawriter.filepattern_fill_slice_number_or_position(filepattern, "")
        base, ext = os.path.splitext(filepattern)
        return base + "_parameters.yaml"

    def __generate_vtk(self, vtk_file="~/tree.vtk"):
        vtk_file = op.expanduser(vtk_file)
        from tree import TreeBuilder

        if "tree_data" in dir(self.gen):
            resolution=self.config["postprocessing"]["measurement_resolution"]
            tvg = TreeBuilder('vtk', generator_params={"cylinder_resolution":resolution, "sphere_resolution":resolution})
            # yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
            # tvg.importFromYaml(yaml_path)
            tvg.voxelsize_mm = self.voxelsize_mm
            tvg.shape = self.gen.areasize_px
            tvg.tree_data = self.gen.tree_data
            output = tvg.buildTree() # noqa
            # tvg.show()
            tvg.saveToFile(vtk_file)
            return tvg.generator.polyData

    def filepattern_fill_series(self):
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
        return filepattern

    def filepattern_split(self):
        filepattern = self.filepattern_fill_series()
        filepattern = re.sub(r"({.*?})", r"", filepattern)
        root, ext = op.splitext(filepattern)
        return root, ext

    def _get_fn_base(self):
        fn_base, fn_ext = self.filepattern_split()

        fn_base = op.expanduser(fn_base)
        return fn_base

    def save_parameters(self):
        if self.parameters_changed_before_save:
            self.run()
        # prepare path to save

        config_filepattern = self.get_config_file_pattern()
        series_number = io3d.datawriter.get_unoccupied_series_number(filepattern=config_filepattern)
        # series_number = io3d.datawriter.get_unoccupied_series_number(filepattern=self.config["filepattern"])
        self.config['filepattern_series_number'] = series_number

        fn_base = self._get_fn_base()

        dirname = op.dirname(fn_base)
        if not op.exists(dirname):
            os.makedirs(dirname)

        io3d.misc.obj_to_file(self.config, filename=fn_base + "_parameters.yaml")

        handler = logging.FileHandler(fn_base + ".log")
        handler.setFormatter(self.formatter)
        handler.setLevel(self.loglevel)
        self.memoryhandler.setTarget(handler)
        self.memoryhandler.flush()
        self.memoryhandler.flushLevel = self.loglevel

    def save_volume(self):
        if self.parameters_changed_before_save:
            self.run()
        # TODO split save_volume and save_parameters
        self.save_parameters()
        import io3d.misc
        import time
        t0 = time.time()
        fn_base = self._get_fn_base()
        # config["filepattern"] = filepattern

        self._numeric_measurement(fn_base)
        self.save_stats(fn_base)
        t1 = time.time()
        logger.debug("before volume generate " + str(t1-t0))
        # postprocessing
        if "generate_volume" in dir(self.gen):
            # self.data3d = self.gen.generate_volume()
            self.data3d = self.gen.generate_volume(dtype="uint8")
            self.voxelsize_mm = self.gen.voxelsize_mm
            postprocessing_params = self.config["postprocessing"]
            data3d = self.postprocessing(**postprocessing_params)
            self.gen.data3d = data3d
        # self.gen.saveVolumeToFile(self.config["filepattern"])
        t2 = time.time()
        logger.debug("before volume save " + str(t2-t0))
        self.gen.saveVolumeToFile(self.filepattern_fill_series())
        t3 = time.time()
        logger.info("time before volume generate: " + str(t1-t0))
        logger.info("time before volume save: " + str(t2-t0))
        logger.info("time after volume save: " + str(t3-t0))
        self.stats_times["numeric_measurement_time_s"] = [t1-t0]
        self.stats_times["generate_volume_time_s"] = [t2-t1]
        self.stats_times["save_volume_time_s"] = [t3-t2]

        # self.memoryhandler.flush()

    def postprocessing(
            self,
            gaussian_blur=True,
            gaussian_filter_sigma_mm=1.0,
            add_noise=True,
            noise_preview=False,
            # gaussian_noise_stddev=10.0,
            # gaussian_noise_center=0.0,
            limit_negative_intensities=True,
            noise_random_generator_seed=0,
            exponent=0.0001,
            lambda_start=0.1,
            lambda_stop=3.0,
            noise_amplitude = 40.0,
            noise_mean = 30.0,
            surface_measurement=False,
            measurement_multiplier=-1,
            measurement_resolution=50,
            output_dtype="uint8"

    ):
        if gaussian_blur:
            sigma_px = gaussian_filter_sigma_mm / self.voxelsize_mm

            self.data3d = scipy.ndimage.filters.gaussian_filter(
                self.data3d,
                sigma=sigma_px)

        if add_noise:

            dt = self.data3d.dtype
            noise = self.generate_noise()
            noise = noise.astype(self.data3d.dtype)
            # noise = np.random.normal(loc=gaussian_noise_center, scale=gaussian_noise_stddev, size=self.data3d.shape)
            self.data3d = (self.data3d + noise).astype(self.data3d.dtype)

        if limit_negative_intensities:
            self.data3d[self.data3d < 0] = 0
        self.config["postprocessing"]["measurement_multiplier"] = measurement_multiplier

        return self.data3d

    def generate_noise(self):
        import ndnoise
        import ndnoise.generator
        pparams = self.config["postprocessing"]
        # data3d = self.postprocessing(**postprocessing_params)
        noise = ndnoise.noises(
            shape=self.gen.areasize_px,
            sample_spacing=self.gen.voxelsize_mm,
            exponent=pparams["exponent"],
            random_generator_seed=pparams["noise_random_generator_seed"],
            lambda_start=pparams["lambda_start"],
            lambda_stop=pparams["lambda_stop"],

        ).astype(np.float16)
        mx = np.max(noise)
        noise = pparams["noise_amplitude"] * noise/mx
        noise += pparams["noise_mean"]
        return noise

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

    def _numeric_measurement(self, fn_base):
        # import numpy as np
        from tree import TreeBuilder
        vxsz = self.config["areasampling"]["voxelsize_mm"]
        vxsz = np.asarray(vxsz).astype(np.float) / self.config["postprocessing"]["measurement_multiplier"]
        shape = self.config["areasampling"]["areasize_px"]
        measurement_multiplier = self.config["postprocessing"]["measurement_multiplier"]
        surface_measurement = self.config["postprocessing"]["surface_measurement"]
        if measurement_multiplier > 0 and surface_measurement:
            shape = np.asarray(shape) * measurement_multiplier
            self._numeric_surface_measurement_shape = shape

            shape = shape.astype(np.int)
            tvgvol = TreeBuilder("vol")
            tvgvol.voxelsize_mm = vxsz
            tvgvol.shape = shape
            tvgvol.tree_data = self.gen.tree_data

            data3d = tvgvol.buildTree()
            import measurement
            surface, vertices, faces = measurement.surface_measurement(data3d, tvgvol.voxelsize_mm, return_vertices_and_faces=True)
            self._surface_vertices = vertices
            self._surface_faces = faces

            volume = np.sum(data3d > 0) * np.prod(vxsz)

            self.dataframes["overall"]["surf. num. [mm^2]"] = [surface]
            self.dataframes["overall"]["vol. num. [mm^3]"] = [volume]

            filename = fn_base + "_raw_{:06d}.jpg"
            import io3d.misc
            data = {
                'data3d': data3d.astype(np.uint8) * 70, #* self.output_intensity,
                'voxelsize_mm': vxsz,
                # 'segmentation': np.zeros_like(self.data3d, dtype=np.int8)
            }
            io3d.write(data, filename)
        else:
            surface = 0
        return surface

    def prepare_stats(self):
        to_rename = {
            "length": "length [mm]",
            "volume": "volume [mm^3]",
            "surface": "surface [mm^2]",
            "radius": "radius [mm]"
        }
        to_rename_density = {
            "length": "length d. [mm^-2]",
            "volume": "volume d. []",
            "surface": "surface d. [mm^-1]"
            # "radius": "radius [mm^-2]"
        }

        df = self.gen.getStats()
        self.dataframes["elements"] = df

        dfdescribe = df.describe()
        dfdescribe.insert(0, "", dfdescribe.index)
        count = dfdescribe["length"][0]
        dfdescribe = dfdescribe.ix[1:]
        dfdescribe = dfdescribe.rename(columns=to_rename)
        self.dataframes["describe"] = dfdescribe

        dfmerne = df[["length", "volume", "surface"]].sum() / self.gen.area_volume
        dfmernef = dfmerne.to_frame().transpose().rename(columns=to_rename_density)
        self.dataframes["density"] = dfmernef

        # whole sumary data
        dfoverall = df[["length", "volume", "surface"]].sum()
        dfoverallf = dfoverall.to_frame().transpose().rename(columns=to_rename)
        dfoverallf["area volume [mm^3]"] = [self.gen.area_volume]
        dfoverallf["count []"] = [count]


        import vtk
        mass = vtk.vtkMassProperties()
        # mass.SetInputData(object1Tri.GetOutput())
        mass.SetInputData(self.polydata)
        surf = mass.GetSurfaceArea()
        vol = mass.GetVolume()
        dfoverallf["numeric volume [mm^2]"]=[vol]
        dfoverallf["numeric surface [mm^2]"]=[surf]
        self.dataframes["overall"] = dfoverallf

    def save_stats(self, fn_base):
        for dfname in self.dataframes:
            self.dataframes[dfname].to_csv(fn_base + "_" + dfname + ".csv")


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
    def run_batch(self, config_list):
        for filename in config_list:
            if filename is None:
                continue
            params = io3d.misc.obj_from_file(filename=filename)
            default_config = self.get_default_config()
            self.update_config(**default_config)
            self.update_config(**params)

            self.run()
            self.save_volume()


class ConfigFileManager():
    def __init__(
            self,
            appname=None,
            config_dir_pattern="~/.config/",
            default_config_file="default_config.yaml",
            favorite_config_file="favorite_config.yaml",
            init_config_file="init_config.yaml",
            log_file="favorite.yaml"
    ):
        self.appname = appname
        self.config_dir = op.expanduser(op.join(config_dir_pattern, appname))

        self.default_config_file = op.join(self.config_dir, default_config_file)
        self.default_config = None
        self.favorite_config_file = op.join(self.config_dir, favorite_config_file)
        self.favorite_config = None
        self.init_config_file = op.join(self.config_dir, init_config_file)
        self.init_config = None
        self.logfile= op.join(self.config_dir, log_file)

    def init_config_dir(self):
        if not op.exists(self.config_dir):
            import os
            os.makedirs(self.config_dir)

    def save_default(self, config):
        io3d.misc.obj_to_file(config, self.default_config_file)

    def load_default(self):
        return io3d.misc.obj_from_file(self.default_config_file)

    def save_favorite(self, config):
        io3d.misc.obj_to_file(config, self.favorite_config_file)

    def load_favorite(self):
        return io3d.misc.obj_from_file(self.favorite_config_file)

    def save_init(self, config):
        io3d.misc.obj_to_file(config, self.init_config_file)

    def load_init(self):
        return io3d.misc.obj_from_file(self.init_config_file)

def main():
    logger = logging.getLogger()

    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    logger.addHandler(ch)

    config_file_manager = ConfigFileManager("teigen")
    config_file_manager.init_config_dir()
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
        default=config_file_manager.init_config_file,
        # required=True,
        help='input parameter file'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='Debug mode')

    parser.add_argument(
        '-ni', '--nointeractivity', action='store_true',
        help='No interactivity mode')

    parser.add_argument(
        '-l', '--logfile',
        default="~/teigen.log",
        help='Debug mode')
    args = parser.parse_args()


    if args.debug:
        ch.setLevel(logging.DEBUG)

    #default param file
    if not op.exists(op.expanduser(args.parameterfile)):
        args.parameterfile = None

    if args.nointeractivity:
        tg = Teigen(logfile=args.logfile)
        if args.parameterfile is not None:
            params = io3d.misc.obj_from_file(args.parameterfile)
            tg.update_config(**params)
        tg.run()
        # tg.run(**params)
        tg.save_volume()
    else:
        app = QApplication(sys.argv)
        params = None
        if args.parameterfile is not None:
            params = io3d.misc.obj_from_file(args.parameterfile)
        cw = TeigenWidget(logfile=args.logfile, config=params)
        cw.show()
        app.exec_()


if __name__ == "__main__":
    main()