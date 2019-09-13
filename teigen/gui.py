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
import PyQt5
from PyQt5.QtWidgets import (QLabel, QPushButton, QGridLayout, QTabWidget,
                         QStatusBar, QFileDialog, QToolBar)

from .tgmain import main
from PyQt5 import QtCore, QtGui, QtWidgets
import os.path as op
import copy
import collections
# from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
# from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import matplotlib.pyplot as plt

from . import dictwidgetqt
from . import iowidgetqt
import io3d.datawriter
import io3d.misc
from . import dictwidgetpg


from pyqtgraph.parametertree import Parameter, ParameterTree

from .tgmain import Teigen, CKEY_OUTPUT, CKEY_APPEARANCE, CKEY_MEASUREMENT
from .teigendoc import teigendoc, teigen_keysdoc


class TeigenWidget(QtWidgets.QWidget):
    def __init__(self, ncols=2, qapp=None, logfile="~/teigen.log", config=None, use_default_config=False):
        super(TeigenWidget, self).__init__()
        self.logfile = logfile
        self.ncols = ncols
        self.gen = None
        self.figures = {}
        self.ui_stats_shown = False
        self.teigen = Teigen(logfile=self.logfile)
        if use_default_config:
            self.teigen.use_default_config()
        if config is not None:
            self.teigen.update_config(**config)
        self.version = self.teigen.version
        self.config = {}
        self.run_number = 0
        self.qapp = qapp
        self.init_ui()

    def collect_config_from_gui(self):
        id = self._ui_generators_tab_wg.currentIndex()
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
        config[CKEY_APPEARANCE] = area_cfg["Appearance"]
        config[CKEY_OUTPUT] = area_cfg["Output"]
        config[CKEY_MEASUREMENT] = area_cfg["Measurement"]
        self.config = config

    def _parameters_changed(self, param, changes):
        logger.debug("parameters changed")
        logger.debug("parameters changed")
        self.teigen.parameters_changed_before_save = True
        self._ui_btn_step2.setEnabled(False)
        # self.on_config_update()

    def collect_config_from_gui_and_push_to_teigen(self):
        self.collect_config_from_gui()
        self.teigen.update_config(**self.config)

    def step1(self):
        self.collect_config_from_gui_and_push_to_teigen()

        # self.config = new_cfg
        self.teigen.step1()

    def _ui_show_potential_output_path(self):
        fn = self.teigen.config["filepattern_abspath"]
        # if fn is None:
        #     fn = self.teigen.filepattern_fill_potential_series()
        self._ui_output_path.setText(fn)
        logger.debug("output path refreshed " + fn)

    def _show_stats_after_step1(self):
        logger.debug("show stats after step1 begin ")
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
            self.mainLayout.addWidget(self.stats_tab_wg, 0, 3, 6, 2)

        self.actual_subtab_wg = QTabWidget()
        self.stats_tab_wg.addTab(self.actual_subtab_wg, '' + run_number_alpha)
        logger.debug("initiating canvas for graphs")
        if True:
            from matplotlib.figure import Figure
            # self.figure = Figure(figsize=(5, 3))
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

        from . import tablewidget
        # import ipdb; ipdb.set_trace()
        # TODO take care about redrawing
        # self.stats_tab_wg.addTab(self._wg_tab_merne, "Density table")
        logger.debug("tabs initiatization")
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

        # TODO move to main column window
        self._wg_btn_tab_save = QPushButton("Save in one row", self)
        self._wg_btn_tab_save.setToolTip("Save all data in one row")
        self._wg_btn_tab_save.clicked.connect(self.btn_save_in_one_row)

        self._wg_tables = QtWidgets.QWidget()
        self._wg_tables.setLayout(QGridLayout())
        self._wg_tables.layout().addWidget(self._wg_tab_describe)
        self._wg_tables.layout().addWidget(self._wg_tab_merne)
        self._wg_tables.layout().addWidget(self._wg_tab_overall)
        self._wg_tables.layout().addWidget(self._wg_btn_tab_save)

        self._wg_tab_describe.show()
        self._wg_tab_describe.raise_()
        # self.mainLayout.addWidget(self._wg_tab_describe, 0, 2, 5, 2)
        # self.stats_tab_wg.addTab(self._wg_tab_describe, "Stats table")
        self.actual_subtab_wg.addTab(self._wg_tables, "Summary " + run_number_alpha)
        # self.resize(600,700)


        logger.debug("poly data visualization init")
        if self.teigen.polydata_volume is not None and self.teigen.config[CKEY_APPEARANCE]["surface_3d_preview"]:
            import imtools.show_segmentation_qt
            logger.debug("segmentation widget loading")
            # test code
            fn = op.expanduser("~/lisa_data/sample_data/liver-seg001.mhd"),
            # datap = io3d.read(fn, dataplus_format=True)
            import imtools.sample_data
            datap = imtools.sample_data.donut()

            segmentation = datap['segmentation']
            voxelsize_mm = datap['voxelsize_mm']
            self._wg_show_3d = imtools.show_segmentation_qt.ShowSegmentationWidget(
                # datap["segmentation"],
                # show_load_button=True
            )
            QtWidgets.QApplication.processEvents()

            logger.debug("read polydata")
            # TODO use again - unstability is not here
            # this ted to be unstable
            # self._wg_show_3d.add_vtk_polydata(self.teigen.polydata_volume)
            # so we are using file way
            temp_vtk_file = op.expanduser(self.teigen.temp_vtk_file)
            self.teigen.save_surface_to_file(temp_vtk_file)
            self._wg_show_3d.add_vtk_file(temp_vtk_file)

            # test codee
            # self._wg_show_3d.add_
            logger.debug("init new tab")
            self.actual_subtab_wg.addTab(self._wg_show_3d, "Visualization " + run_number_alpha)
        else:
            self._wg_show_3d = None

        self.ui_stats_shown = True

        logger.debug("noise preview init")
        if (
                    self.teigen.config[CKEY_APPEARANCE]["noise_preview"] and
                    self.teigen.config["postprocessing"]["add_noise"]):
            self._noise_figure = plt.figure()
            self._noise_canvas = FigureCanvas(self._noise_figure)
            # self.toolbar = NavigationToolbar(self.canvas, self)
            self.actual_subtab_wg.addTab(self._noise_canvas, 'Noise ' + run_number_alpha)
            noise = self.teigen.generate_noise()
            plt.imshow(noise[0, :, :], cmap="gray")
            plt.colorbar()

        logger.debug("show potential output path")
        self._ui_show_potential_output_path()

    def update_stats(self):
        """
        Function is used after volume generation and save.
        :return:
        """
        from . import tablewidget
        self._wg_tab_overall.deleteLater()
        self._wg_tab_overall = None
        #     self._wg_tab_describe.deleteLater()
        # self._wg_tables
        dfoverall = self.teigen.dataframes["overall"]
        self._wg_tab_overall = tablewidget.TableWidget(self, dataframe=dfoverall)
        self._wg_tab_overall.setMaximumHeight(80)
        self._wg_tables.layout().addWidget(self._wg_tab_overall)

        output_path = QLabel()
        output_path.setText(self.teigen.filepattern_fill_series())
        self._wg_tables.layout().addWidget(output_path)
        #     self._wg_tab_describe.deleteLater()
        #     self._wg_tab_describe = None
        #     self._wg_tab_merne.deleteLater()
        #     self._wg_tab_merne = None

        # Show surface
        measurement_multiplier = self.teigen.config[CKEY_OUTPUT]["aposteriori_measurement_multiplier"]
        surface_measurement = self.teigen.config[CKEY_OUTPUT]["aposteriori_measurement"]
        show_surface = self.teigen.config[CKEY_APPEARANCE]["show_aposteriori_surface"]
        if surface_measurement and (measurement_multiplier > 0) and show_surface:
            fig = plt.figure()
            self._surface_figure = fig
            self._surface_canvas = FigureCanvas(self._surface_figure)
            # self.toolbar = NavigationToolbar(self.canvas, self)
            run_number_alpha = chr(ord("A") + self.run_number)
            self.actual_subtab_wg.addTab(self._surface_canvas, 'Aposteriori Surface ' + run_number_alpha)

            # import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
            faces, vertices = self.teigen.get_aposteriori_faces_and_vertices()
            # fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111, projection='3d')

            # Fancy indexing: `verts[faces]` to generate a collection of triangles
            mesh = Poly3DCollection(vertices[faces])
            mesh.set_edgecolor('r')
            ax.add_collection3d(mesh)
            sh = self.teigen._numeric_surface_measurement_shape
            ax.set_xlim(0, sh[0])  # a = 6 (times two for 2nd ellipsoid)
            ax.set_ylim(0, sh[1])  # b = 10
            ax.set_zlim(0, sh[2])  # c = 16
            # plt.show()

    def complicated_to_yaml(self, cfg):
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

    def _get_generator(self, id):
        pass

    def init_ui(self):
        wtitle = "Teigen " + self.version
        self.setWindowTitle(wtitle)
        self.mainLayout = QGridLayout(self)

        self.statusBar = QStatusBar()
        self.mainLayout.addWidget(self.statusBar, 10, 0, 1, 2)
        self.progressBar = QtWidgets.QProgressBar()
        self.progressBar.setRange(0, 10000)
        self.progressBar.setValue(0)

        self.ui_stop_button = QPushButton("Stop", self)
        self.ui_stop_button.clicked.connect(self.btnStop)

        self.statusBar.addWidget(self.progressBar)
        self.statusBar.addWidget(self.ui_stop_button)
        self.progressBar.show()

        self.configBarLayout = QGridLayout(self)
        self._ui_config_init()


    def _ui_init_buttons(self):

        # Toolbar
        # Use default config
        self._ui_btn_default_config = QtWidgets.QAction(self.style().standardIcon(QtWidgets.QStyle.SP_BrowserReload), "Reset parameters to default values", self)
        self._ui_btn_default_config.setToolTip("Reset all parameters to default value (Ctrl+R)")
        self._ui_btn_default_config.setShortcut('Ctrl+R')
        self._ui_btn_default_config.triggered.connect(self.btn_use_default_config)

        # Load
        self._ui_btn_load_config = QtWidgets.QAction(self.style().standardIcon(QtWidgets.QStyle.SP_DialogOpenButton), "Load params", self)
        self._ui_btn_load_config.setToolTip("Load params from file with file dialog (Ctrl+L)")
        self._ui_btn_load_config.setShortcut('Ctrl+L')
        self._ui_btn_load_config.triggered.connect(self.btn_load_config)
        # self.configBarLayout.addWidget(self._ui_btn_load_config, 1, 3, 1, 1)  # , (gd_max_i / 2), text_col)

        self._ui_btn_save = QtWidgets.QAction(self.style().standardIcon(QtWidgets.QStyle.SP_DialogSaveButton), "Save parameters", self)
        self._ui_btn_save.setToolTip("Save generator parameters (Ctrl+S)")
        self._ui_btn_save.triggered.connect(self.btn_save_parameters)
        self._ui_btn_save.setShortcut('Ctrl+S')
        # self.configBarLayout.addWidget(self._ui_btn_save, 1, 1, 1, 1)

        self._ui_btn_save_and_add_to_batch = QtWidgets.QAction(self.style().standardIcon(QtWidgets.QStyle.SP_FileDialogDetailedView), "Save parameters and add to batch", self)
        self._ui_btn_save_and_add_to_batch.setToolTip("Save generator parameters and then add to batch (Ctrl+B)")
        self._ui_btn_save_and_add_to_batch.triggered.connect(self.btn_save_parameters_and_add_to_batch)
        self._ui_btn_save_and_add_to_batch.setShortcut('Ctrl+B')
        # self.configBarLayout.addWidget(self._ui_btn_save_and_add_to_batch, 1, 2, 1, 1)  # , (gd_max_i / 2), text_col)

        # exitAction = QtGui.QAction(QtGui.QIcon('exit24.png'), 'Exit', self)
        # self.exitAction = QtGui.QAction( self.style().standardIcon(QtGui.QStyle.SP_DialogOpenButton) , 'Exit', self)
        # self.exitAction.setShortcut('Ctrl+Q')
        # self.exitAction.triggered.connect(QtGui.qApp.quit)

        # loadTubeSkeletonAction = QtGui.QAction( QtCore.QCoreApplication.translate("ahoj"), 'Exit', self)
        self.loadTubeSkeletonAction = QtWidgets.QAction(self.style().standardIcon(QtWidgets.QStyle.SP_MediaSkipForward), "Load skeleton",  self)
        self.loadTubeSkeletonAction.setShortcut('Ctrl+K')
        self.loadTubeSkeletonAction.setToolTip('Make Step 1 by loading skeleton (Ctrl+K)')
        self.loadTubeSkeletonAction.triggered.connect(self.btn_step1_by_load_tube_skeleton)

        self.toolbar = QToolBar("File")
        self.mainLayout.addWidget(self.toolbar, 3, 1, 1, 1)
        self.toolbar.addAction(self._ui_btn_save)
        self.toolbar.addAction(self._ui_btn_save_and_add_to_batch)
        self.toolbar.addAction(self._ui_btn_default_config)
        self.toolbar.addAction(self._ui_btn_load_config)
        self.toolbar.addAction(self.loadTubeSkeletonAction)
        # self.toolbar.addAction(self.exitAction)
        self._ui_output_path = QLabel(self)
        self._ui_output_path.setText("")



        self.mainLayout.addWidget(self._ui_output_path, 2, 1, 1, 2)  # , (gd_max_i / 2), text_col)
        # self._ui_show_potential_output_path()

        self.mainLayout.addLayout(self.configBarLayout, 3, 1, 1, 2)  # , (gd_max_i / 2), text_col)

        self._ui_btn_step1 = QPushButton("Step 1 - Preview - Generate skeleton", self)
        self._ui_btn_step1.setToolTip("Generate preview and skeleton (Ctrl+P)")
        self._ui_btn_step1.clicked.connect(self.btnRunStep1)
        self._ui_btn_step1.setShortcut('Ctrl+P')
        self.mainLayout.addWidget(self._ui_btn_step1, 4, 1, 1, 2)  # , (gd_max_i / 2), text_col)

        # self.posprocessing_wg = dictwidgetqt.DictWidget(postprocessing_params)
        # self.mainLayout.addWidget(self.posprocessing_wg, 3, 1)

        self._ui_btn_step2 = QPushButton("Step 2 - Generate and save volumetric data", self)
        self._ui_btn_step2.setToolTip("Save image slices and meta information (Ctrl+R)")
        self._ui_btn_step2.clicked.connect(self.btnRunStep2)
        self._ui_btn_step2.setShortcut('Ctrl+R')
        self.mainLayout.addWidget(self._ui_btn_step2, 5, 1, 1, 2)  # , (gd_max_i / 2), text_col)

    def _ui_config_init(self):
        self._ui_init_buttons()
        self.ui_output_dir_widget = iowidgetqt.SetDirWidget(
            self.teigen.config["filepattern"], "output directory")
        self.ui_output_dir_widget.setToolTip(
            "Data are stored in defined directory.\nOutput format is based on file extension.\n" +
            "For saving into image stack use 'filename{:06d}.jpg'")
        self.mainLayout.addWidget(self.ui_output_dir_widget, 1, 1, 1, 2)  # , (gd_max_i / 2), text_col)

        postprocessing_params = self.teigen.config["postprocessing"]

        hide_keys = ["build", "gtree", "voxelsize_mm", "areasize_px",
                     "resolution", "n_slice", "dims", "tube_shape",
                     "radius_distribution_normal", "radius_distribution_uniform",
                     "radius_distribution_fixed", "allow_overlap"
                     ]
        self._ui_generators_tab_wg = QTabWidget()
        self._ui_generators_tab_wg.setMinimumWidth(400)
        self.mainLayout.addWidget(self._ui_generators_tab_wg, 0, 1, 1, 2)
        # l = QVBoxLayout(self)

        rename_captions_dict = {
            "voxelsize_mm": "voxel size [mm]",
        }
        dropdownoboxes = {
            "radius_distribution": ["normal", "fixed", "uniform"]
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
                dropdownboxes=dropdownoboxes,
            )
            self._ui_generator_widgets.append(wg)
            self._ui_generators_tab_wg.addTab(wg, generator_name.replace(" ", "\n"))
            # self._ui_generators_tab_wg.

        # l.addWidget(self._ui_generators_tab_wg)
        # wgn = QtGui.QWidget()
        # layout = QtGui.QFormLayout()
        # layout.addRow("Name", QtGui.QLineEdit())
        # layout.addRow("Address",QtGui.QLineEdit())
        # wgn.setLayout(layout)
        # self._ui_generators_tab_wg.addTab(wgn, "ahoj")
        id = self.teigen.get_generator_id_by_name_or_number(teigen_config["generator_id"])
        self._ui_generators_tab_wg.setCurrentIndex(id)
        # self.mainLayout.setColumnMinimumWidth(text_col, 500)

        # pyqtgraph experiments
        input_params = {
            "Area Sampling": dictwidgetpg.AreaSamplingParameter(name='Area Sampling',
                                                                **self.teigen.config["areasampling"]),
            "Postprocessing": postprocessing_params,
            "Batch processing": dictwidgetpg.BatchFileProcessingParameter(
                name="Batch processing", children=[
                    {'name': 'Run batch', 'type': 'action', "tip": "output"},
                ]),
            "Appearance": self.teigen.config["appearance"],
            "Output": self.teigen.config["output"],
            "Measurement": self.teigen.config["measurement"]
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
        dictwidgetpg.add_tip(gr_struct, "noise_preview", "this is noise")
        dictwidgetpg.add_tips(gr_struct, teigen_keysdoc)
        gr_struct["children"][1]['tip'] = "Apperance tip"
        # import ipdb; ipdb.set_trace()
        # gr_struct["children"][2]['tip'] = "output tip"
        gr_struct["children"][3]['tip'] = "post processing tip"
        gr_struct["children"][4]['tip'] = "measurement"
        p = Parameter.create(**gr_struct)

        t = ParameterTree()
        t.setParameters(p, showTop=False)
        t.setMinimumWidth(380)
        # t.setColumnCount(3)
        t.show()

        p.sigTreeStateChanged.connect(self._parameters_changed)
        p.param('Batch processing', 'Run batch').sigActivated.connect(self.run_batch)

        # how to add button
        # i5  = pg.TreeWidgetItem(["Item 5"])
        # b5 = QtGui.QPushButton('Button')
        # i5.setWidget(1, b5)
        # t.addTopLevelItem(i5)

        self.mainLayout.addWidget(t, 0, 0, 6, 1)
        self.config_wg = t
        # self.config_wg.setToolTip(teigendoc)
        self.area_sampling_params = p
        self.teigen.progress_callback = self._progressbar_update
        self._ui_btn_step2.setEnabled(False)

    def delete_wg(self, wg):
        self.mainLayout.removeWidget(wg)
        wg.deleteLater()
        wg = None

    def _ui_config_deinit(self):

        self.delete_wg(self.ui_output_dir_widget)
        self.delete_wg(self.config_wg)
        self.delete_wg(self._ui_generators_tab_wg)
        self.delete_wg(self._ui_btn_step1)
        self.delete_wg(self._ui_btn_step2)
        self.delete_wg(self.toolbar)
        # self.delete_wg(self._ui_btn_save_and_add_to_batch)
        # self.delete_wg(self._ui_btn_save)

    def step1_by_load_tube_skeleton(self, filename):
        self.collect_config_from_gui_and_push_to_teigen()
        self.teigen.step1_by_load_tube_skeleton(filename)
        self._show_stats_after_step1()
        self._ui_btn_step2.setEnabled(True)

    def btn_step1_by_load_tube_skeleton(self):
        fn = op.dirname(self.teigen.get_fn_base())
        filename = QFileDialog.getOpenFileName(self, 'Open config file',
                                               fn, "Config files (*.yaml)")[0]
        if filename is not None:
            filename = str(filename)
        self.step1_by_load_tube_skeleton(filename=filename)

    def btn_load_config(self):
        self.load_config()

    def btn_use_default_config(self):
        self.teigen.use_default_config()
        self._ui_config_deinit()
        self._ui_config_init()

    def load_config(self, filename=None):
        """
        load config from file, if filename is none, file dialog is created
        :param filename:
        :return:
        """
        if filename is None:
            fn = op.dirname(self.teigen.get_fn_base())
            filename = QFileDialog.getOpenFileName(self, 'Open config file',
                                                   fn, "Config files (*.yaml)")[0]
            if filename is not None:
                filename = str(filename)
        if filename is not None:
            params = io3d.misc.obj_from_file(filename)
            self.teigen.update_config(**params)
            self._ui_config_deinit()
            self._ui_config_init()

    def _progressbar_update(self, obj, level, *args, **kwargs):
        self.progressBar.setValue(int(10000 * level))
        if "statusbar_text" in kwargs:
            # add this in gui
            logger.debug("statusbar_text " + kwargs["statusbar_text"])
            ## end of pyqtgraph tree

    def run_batch(self):

        none, area_cfg = dictwidgetpg.from_pyqtgraph_struct(self.area_sampling_params.saveState())
        lst = area_cfg["Batch processing"].values()
        self.teigen.run_batch(lst)

    def btn_save_parameters(self):
        self.save_parameters()

    def save_parameters(self, filename=None):
        if filename is None:
            # fn = op.dirname(self.teigen.get_fn_base())
            fn = op.dirname(self.teigen.get_fn_base())
            fn += "_parameters.yaml"
            filename = QFileDialog.getSaveFileName(self, 'Save config file',
                                                   fn, "Config files (*.yaml)")[0]
            if filename is not None:
                filename = str(filename)
        self.collect_config_from_gui_and_push_to_teigen()
        self.teigen.save_parameters(filename=filename)
        return filename

    def btn_save_parameters_and_add_to_batch(self):
        self.save_parameters_and_add_to_batch()

    def save_parameters_and_add_to_batch(self, filename=None):
        self.collect_config_from_gui_and_push_to_teigen()
        config_filename = self.teigen.save_parameters(filename=filename)
        self.area_sampling_params.param("Batch processing").add_filename(config_filename)

    def btnRunStep1(self):

        logger.debug("btnRunStage1")
        # logger.debug(str(self.config))
        self.step1()
        self._show_stats_after_step1()
        self.run_number += 1
        self._ui_btn_step2.setEnabled(True)

    def btnStop(self):
        self.teigen.stop()
        pass

    def btn_save_in_one_row(self):
        # fn = op.dirname(self.teigen.get_fn_base())
        fn_base, fn_ext = self.teigen.filepattern_split()
        fn = op.dirname(fn_base)
        superior_dir, filen = op.split(fn)
        fn = op.join(superior_dir, "output_rows.csv")
        filename = QFileDialog.getSaveFileName(self,
                                               'Save config file',
                                               fn,
                                               "(*.csv)",
                                               options=QtWidgets.QFileDialog.DontConfirmOverwrite
                                               )[0]
        text, ok = QtWidgets.QInputDialog.getText(self, 'Note Dialog',
                                              'Note:')
        if filename is not None:
            filename = str(filename)
            self.teigen.save_stats_to_row(filename, note=text)

    def btnRunStep2(self):
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
            self.step1()
            self._show_stats_after_step1()

        # filename = op.join(self.ui_output_dir_widget.get_dir(), filename)
        filename = self.config["filepattern"]

        # filename = iowidgetqt.str_format_old_to_new(filename)

        self.teigen.step2()
        fn_base, fn_ext = self.teigen.filepattern_split()
        fn_base = self.teigen.config["filepattern_abspath"]
        self.figure.savefig(fn_base + "_" + "graph.pdf")
        self.figure.savefig(fn_base + "_" + "graph.png")
        self.figure.savefig(fn_base + "_" + "graph.svg")
        self.figure.savefig(fn_base + "_" + "graph.eps")

        from PyQt5.QtGui import QPixmap
        # from PyQt5.QtGui import Q

        if self._wg_show_3d is not None:
            # p = QPixmap.grabWidget(self._wg_show_3d.vtkWidget)
            p = self._wg_show_3d.vtkWidget.grab()
            p.save(fn_base + "_snapshot.png", 'png')

        # self.teigen.gen.saveVolumeToFile(filename)
        self.update_stats()

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


if __name__ == "__main__":
    main()
