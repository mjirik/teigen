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
# conda install -c conda-forge begins
# import begin

import PyQt4
from PyQt4.QtGui import QLabel, \
    QPushButton, QGridLayout, QTabWidget, QStatusBar, \
    QFileDialog

from tgmain import main

from PyQt4 import QtGui
import os.path as op
import copy
import collections

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
# from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import matplotlib.pyplot as plt

import dictwidgetqt
import iowidgetqt
import io3d.datawriter
import io3d.misc
import dictwidgetpg
from pyqtgraph.parametertree import Parameter, ParameterTree

from tgmain import CKEY_OUTPUT, CKEY_APPEARANCE, Teigen


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

    def _ui_show_potential_output_path(self):
        fn = self.teigen.filepattern_fill_potential_series()
        self._ui_output_path.setText(fn)
        logger.debug("output path refreshed " + fn)

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

        # TODO move to main column window
        self._wg_btn_tab_save = QPushButton("Save in one row", self)
        self._wg_btn_tab_save.setToolTip("Save all data in one row")
        self._wg_btn_tab_save.clicked.connect(self.btn_save_in_one_row)

        self._wg_tables = QtGui.QWidget()
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


        if self.teigen.polydata is not None:
            import imtools.show_segmentation_qt
            self._wg_show_3d = imtools.show_segmentation_qt.ShowSegmentationWidget(None, show_load_button=False)

            # self._wg_show_3d.add_vtk_file(op.expanduser(self.teigen.temp_vtk_file))
            self._wg_show_3d.add_vtk_polydata(self.teigen.polydata)
            self.actual_subtab_wg.addTab(self._wg_show_3d, "Visualization " + run_number_alpha)

        self.ui_stats_shown = True

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

        self._ui_show_potential_output_path()

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
        if surface_measurement and measurement_multiplier > 0 and show_surface:
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

        self.configBarLayout = QGridLayout(self)

        self._ui_btn_load_config = QPushButton("Load params", self)
        self._ui_btn_load_config.setToolTip("Load params from file with file dialog")
        self._ui_btn_load_config.clicked.connect(self.btn_load_config)
        self.configBarLayout.addWidget(self._ui_btn_load_config, 1, 3, 1, 1)  # , (gd_max_i / 2), text_col)

        self._ui_btn_save = QPushButton("Save parameters", self)
        self._ui_btn_save.setToolTip("Save generator parameters")
        self._ui_btn_save.clicked.connect(self.btn_save_parameters)
        self.configBarLayout.addWidget(self._ui_btn_save, 1, 1, 1, 1)

        self._ui_btn_save_and_add_to_batch = QPushButton("Save parameters and add to batch", self)
        self._ui_btn_save_and_add_to_batch.setToolTip("Save generator parameters and then add to batch")
        self._ui_btn_save_and_add_to_batch.clicked.connect(self.btn_save_parameters_and_add_to_batch)
        self.configBarLayout.addWidget(self._ui_btn_save_and_add_to_batch, 1, 2, 1, 1)  # , (gd_max_i / 2), text_col)

        self._ui_output_path = QLabel(self)
        self._ui_output_path.setText("")
        self.mainLayout.addWidget(self._ui_output_path, 2, 1, 1, 2)  # , (gd_max_i / 2), text_col)
        self._ui_show_potential_output_path()

        self.mainLayout.addLayout(self.configBarLayout, 3, 1, 1, 2)  # , (gd_max_i / 2), text_col)

        self._ui_btn_step1 = QPushButton("Step 1 - Preview - Generate skeleton", self)
        self._ui_btn_step1.clicked.connect(self.btnRunStep1)
        self.mainLayout.addWidget(self._ui_btn_step1, 4, 1, 1, 2)  # , (gd_max_i / 2), text_col)

        # self.posprocessing_wg = dictwidgetqt.DictWidget(postprocessing_params)
        # self.mainLayout.addWidget(self.posprocessing_wg, 3, 1)

        self._ui_btn_step2 = QPushButton("Step 2 - Generate and save volumetric data", self)
        self._ui_btn_step2.setToolTip("Save image slices and meta information")
        self._ui_btn_step2.clicked.connect(self.btnRunStep2)
        self.mainLayout.addWidget(self._ui_btn_step2, 5, 1, 1, 2)  # , (gd_max_i / 2), text_col)

        self._ui_config_init()

    def _ui_config_init(self):

        self.ui_output_dir_widget = iowidgetqt.SetDirWidget(
            self.teigen.config["filepattern"], "output directory")
        self.ui_output_dir_widget.setToolTip(
            "Data are stored in defined directory.\nOutput format is based on file extension.\n\
For saving into image stack use 'filename{:06d}.jpg'")
        self.mainLayout.addWidget(self.ui_output_dir_widget, 1, 1, 1, 2)  # , (gd_max_i / 2), text_col)

        postprocessing_params = self.teigen.config["postprocessing"]

        hide_keys = ["build", "gtree", "voxelsize_mm", "areasize_px", "resolution", "n_slice", "dims"]
        self._ui_generators_tab_wg = QTabWidget()
        self.mainLayout.addWidget(self._ui_generators_tab_wg, 0, 1, 1, 2)

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
            self._ui_generators_tab_wg.addTab(wg, generator_name)
        self._ui_generators_tab_wg.setCurrentIndex(teigen_config["generator_id"])
        # self.mainLayout.setColumnMinimumWidth(text_col, 500)

        # pyqtgraph experiments
        input_params = {
            "Area Sampling": dictwidgetpg.AreaSamplingParameter(name='Area Sampling',
                                                                **self.teigen.config["areasampling"]),
            "Postprocessing": postprocessing_params,
            "Batch processing": dictwidgetpg.BatchFileProcessingParameter(
                name="Batch processing", children=[
                    {'name': 'Run batch', 'type': 'action'},
                ]),
            "Appearance": self.teigen.config["appearance"],
            "Output": self.teigen.config["output"]
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
        self.config_wg = t
        self.area_sampling_params = p
        self.teigen.progress_callback = self._progressbar_update

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
        self.delete_wg(self._ui_btn_save)
        self.delete_wg(self._ui_btn_save_and_add_to_batch)

    def btn_load_config(self):
        self.load_config()

    def load_config(self, filename=None):
        """
        load config from file, if filename is none, file dialog is created
        :param filename:
        :return:
        """
        if filename is None:
            fn = op.dirname(self.teigen.get_fn_base())
            filename = QFileDialog.getOpenFileName(self, 'Open config file',
                                                   fn, "Config files (*.yaml)")
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
            print "statusbar_text " + kwargs["statusbar_text"]
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
                                                   fn, "Config files (*.yaml)")
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
        self.run()
        self._show_stats()
        self.run_number += 1

    def btnStop(self):
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
                                               options=QtGui.QFileDialog.DontConfirmOverwrite
                                               )
        text, ok = QtGui.QInputDialog.getText(self, 'Note Dialog',
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
        self.figure.savefig(fn_base + "_" + "graph.eps")

        from PyQt4.QtGui import QPixmap
        p = QPixmap.grabWidget(self._wg_show_3d.vtkWidget)
        p.save(fn_base + "_snapshot.png", 'png')

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


if __name__ == "__main__":
    main()
