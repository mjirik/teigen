#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © %YEAR%  <>
#
# Distributed under terms of the %LICENSE% license.

import logging
logger = logging.getLogger(__name__)
import logging.handlers
import argparse
# import begin
import sys
import os
import os.path as op
import inspect
import numpy as np
import scipy
import re
import datetime
import copy
import collections
import pandas as pd
# from . import generators
from .generators import cylinders
from .generators import gensei_wrapper
# import .generators.cylinders
# from .generators.gensei_wrapper
from .generators import unconnected_cylinders
from imtools.dili import get_default_args
from imma import dili
import io3d.datawriter
import io3d.misc
import ndnoise
import ndnoise.generator
# from . import dictwidgetqt
# from . import geometry3d as g3


CKEY_APPEARANCE = "appearance"
CKEY_OUTPUT = "output"
CKEY_MEASUREMENT = "measurement"


class Teigen():
    def __init__(self, logfile='~/tegen.log', loglevel=logging.DEBUG):
        self.config_file_manager = ConfigFileManager("teigen")
        self.config_file_manager.init_config_dir()
        # self.loglevel = loglevel

        self.logger = logging.getLogger()
        logging.basicConfig()
        self.filehandler = logging.handlers.RotatingFileHandler(
            op.expanduser(logfile),
            maxBytes=1000000,
            backupCount=9
        )
        # self.filehandler.setLevel(self.loglevel)
        # formatter = logging.Formatter('%(asctime)s %(name)-18s %(levelname)-8s %(message)s')
        self.formatter = logging.Formatter(
            '%(asctime)s %(levelname)-8s %(name)-18s %(lineno)-5d %(funcName)-12s %(message)s')
        self.filehandler.setFormatter(self.formatter)
        logger.addHandler(self.filehandler)

        # self.memoryhandler = logging.handlers.MemoryHandler(1024*10, logging.DEBUG, streamhandler)
        self.memoryhandler = logging.handlers.MemoryHandler(1024 * 100)  # , logging.DEBUG, streamhandler)
        # self.memoryhandler.setLevel(self.loglevel)

        self.streamhandler = logging.StreamHandler()
        self.streamhandler.setFormatter(self.formatter)

        self.set_loglevel(loglevel)

        logger.info("Starting Teigen")

        self.logfile = logfile
        self.version = "0.3.0"
        self.data3d = None
        self.voxelsize_mm = None
        self.need_run = True
        self.gen = None
        self.generators_classes = [
            # generators.cylinders.CylinderGenerator,
            # generators.gensei_wrapper.GenseiGenerator,
            # generators.cylinders.CylinderGenerator,
            unconnected_cylinders.UnconnectedCylinderGenerator,
            unconnected_cylinders.UnconnectedCylinderGenerator,
            unconnected_cylinders.UnconnectedCylinderGenerator,
            unconnected_cylinders.UnconnectedCylinderGenerator,
        ]
        self.generators_names = [
            # "Voronoi tubes",
            # "Gensei",
            # "Continuous tubes",
            "Unconnected tubes",
            "Connected tubes",
            "Unconnected porosity",
            "Connected porosity",
        ]
        self._cfg_export_fcn = [
            # self._config2generator_general_export,
            # self._config2generator_gensei_export,
            # self._config2generator_general_export,
            self._config2generator_tubes_export,
            self._config2generator_tubes_export,
            self._config2generator_porosity_export,
            self._config2generator_porosity_export,
        ]
        self._cfg_negative = [
            False, False, True, True
        ]
        self.use_default_config()
        self.progress_callback = None
        self.temp_vtk_file = op.expanduser("~/tree.vtk")
        # 3D visualization data, works for some generators
        self.polydata_volume = None
        self.dataframes = {}
        self.stats_times = {
            "datetime": str(datetime.datetime.now())
        }
        self.parameters_changed_before_save = True
        self.fig_3d_render_snapshot = None
        self.tube_skeleton = {}

    def __del__(self):
        self.filehandler.close()

    def set_loglevel(self, loglevel):
        self.loglevel = loglevel
        self.logger.setLevel(self.loglevel)
        self.filehandler.setLevel(self.loglevel)
        self.memoryhandler.setLevel(self.loglevel)
        self.streamhandler.setLevel(self.loglevel)

    def use_default_config(self):
        self.config = self.get_default_config()

    def get_default_config(self):
        """
        Create default configuration.

        Configuration is composed from
        :return:
        """
        config = collections.OrderedDict()
        # self.config["generators"] = [dictwidgetqt.get_default_args(conf) for conf in self.generators_classes]
        hide_keys = ["build", "gtree", "voxelsize_mm", "areasize_px", "resolution",
                     "n_slice", "dims", "intensity_profile_intensity", "intensity_profile_radius"]
        config["generators"] = collections.OrderedDict()
        for generator_cl, generator_name in zip(
                self.generators_classes,
                self.generators_names
        ):
            generator_params = get_default_args(generator_cl)
            generator_params = dili.kick_from_dict(generator_params, hide_keys)
            config["generators"][generator_name] = generator_params

        config["generators"]["Unconnected tubes"]["allow_overlap"] = False
        config["generators"]["Connected tubes"]["allow_overlap"] = True
        config["generators"]["Unconnected porosity"]["allow_overlap"] = False
        config["generators"]["Connected porosity"]["allow_overlap"] = True

        # self.config["generator_id"] = self.generators_names[0]
        config["generator_id"] = 0
        # self.config = self.configs[0]
        config["postprocessing"] = get_default_args(self.postprocessing)
        config["postprocessing"]["intensity_profile_radius"] = [0.4, 0.7, 1.0, 1.3]
        config["postprocessing"]["intensity_profile_intensity"] = [195, 190, 200, 30]
        # config["postprocessing"][""] = dictwidgetqt.get_default_args(self.postprocessing)
        config["areasampling"] = {
            "voxelsize_mm": [0.01, 0.01, 0.01],
            "areasize_mm": [2.0, 2.0, 2.0],
            "areasize_px": [200, 200, 200]
        }
        config["filepattern"] = "~/teigen_data/{seriesn:03d}/data{:06d}.jpg"
        config["filepattern_abspath"] = None
        # config['filepattern_series_number'] = series_number
        # self.config["voxelsize_mm"] = [1.0, 1.0, 1.0]
        # self.config["areasize_mm"] = [100.0, 100.0, 100.0]
        # self.config["areasize_px"] = [100, 100, 100]
        config[CKEY_APPEARANCE] = {
            "show_aposteriori_surface": True,
            "skip_volume_generation": False,
            "noise_preview": False,
            "surface_3d_preview": False,
            "force_rewrite": False,
        }
        # "force_rewrite"  if series number is used on output dir
        config["output"] = {
            "one_row_filename": "~/teigen_data/output_rows.csv",
            "aposteriori_measurement": False,
            "aposteriori_measurement_multiplier": 1.0,
            "note": ""
        }

        config["measurement"] = {
            "polygon_radius_selection_method": "best",
            # "polygon_radius_selection_method": "inscribed"
            "tube_shape": True,
        }
        return config

    def update_config(self, **config):

        if "required_teigen_version" in config.keys():
            reqired_version = config["required_teigen_version"]
            if reqired_version != self.version:
                logger.error(
                    "Wrong teigen version. Required: " + reqired_version + " , actual " + self.version)
                return
        config = copy.deepcopy(config)
        # there can be stored more than our config. F.e. some GUI dict reconstruction information
        self.config = dili.recursive_update(self.config, config)
        self.voxelsize_mm = np.asarray(self.config["areasampling"]["voxelsize_mm"])
        self.areasize_px = np.asarray(self.config["areasampling"]["areasize_px"])
        self.parameters_changed_before_save = True

    def get_generator_id_by_name_or_number(self, id):

        # if id is not nuber but name of generator
        if type(id) == str:
            for i in range(len(self.generators_names)):
                if id == self.generators_names[i]:
                    id = i
                    break

        if type(id) == str:
            logger.error("Unknown generator name: " + id)

        return id




    def _step1_init_generator(self, tube_skeleton=None):
        t0 = datetime.datetime.now()
        logger.info("step1_init_datetime" + str(t0))
        st0 = str(t0)
        self.stats_times["step1_init_datetime"] = st0
        config = copy.deepcopy(self.config)
        self.config_file_manager.save_init(self.config)

        id = config.pop('generator_id')
        id = self.get_generator_id_by_name_or_number(id)
        self.stop_flag = False

        # area_dct = config["areasampling"]
        # area_cfg = self._cfg_export_fcn[id](area_dct)

        area_cfg = self._cfg_export_fcn[id](config)

        # TODO probably unused
        config.update(area_cfg)

        generator_class = self.generators_classes[id]
        # self.config = get_default_args(generator_class)

        # select only parameters for generator
        # generator_default_config = dictwidgetqt.get_default_args(generator_class)
        # generator_config = dictwidgetqt.subdict(config["generators"][id], generator_default_config.keys())
        generator_config = list(config["generators"].items())[id][1]
        generator_config.update(area_cfg)
        # import ipdb;ipdb.set_trace()
        self.gen = generator_class(**generator_config)
        if id == 2:
            self.gen.MAKE_IT_SHORTER_CONSTANT = 0.0
            self.gen.OVERLAPS_ALOWED = True
        self.gen.progress_callback = self.progress_callback
        if tube_skeleton is not None:
            self.gen.tree_data = tube_skeleton
        logger.debug("step1 init generator finished")
        logger.debug("step1 init generator finished")
        return t0

    def _step1_deinit_save_stats(self, t0):
        self.tube_skeleton = self.gen.tree_data
        # import ipdb;ipdb.set_trace()

        t1 = datetime.datetime.now()
        logger.debug("1D structure is generated")
        logger.debug("before vtk generation")
        pdatas = self.__generate_vtk(self.temp_vtk_file)
        logger.debug("generate vtk finished")
        self.polydata_volume = pdatas[0]
        self.polydata_surface = pdatas[1]
        t2 = datetime.datetime.now()
        self.stats_times["step1_total_time_s"] = (t2 - t0).total_seconds()
        self.stats_times["step1_generate_time_s"] = (t1 - t0).total_seconds()
        self.stats_times["step1_generate_vtk_time_s"] = (t2 - t1).total_seconds()
        self.stats_times["step1_finished"] = True
        self.stats_times["step2_finished"] = False
        self.time_run = t2 - t0
        # self.prepare_stats()
        self.config["filepattern_abspath"] = self.filepattern_fill_series()
        one_row_filename = self.config["output"]["one_row_filename"]
        if one_row_filename != "":
            # self.prepare_stats()
            self.save_stats_to_row(one_row_filename)
        else:
            self.prepare_stats()
        logger.info("time: " + str(self.time_run))
        self.need_run = False
        self.parameters_changed_before_save = False

    def step1_by_load_tube_skeleton(self, filename):
        logger.debug("step1_by_loda_tube_skeleton")
        self.load_tube_skeleton(filename=filename)
        logger.debug("tube skeleton loaded")
        t0 = self._step1_init_generator(self.tube_skeleton)
        logger.debug("generator initiated")
        logger.debug("generator initiated")
        # t0 = datetime.datetime.now()
        # st0 = str(t0)
        # logger.info("step1_init_datetime " + st0)
        # self.stats_times["step1_init_datetime"] = st0
        self._step1_deinit_save_stats(t0)

    def step1(self):
        t0 = self._step1_init_generator()

        # self.gen = generators.gensei_wrapper.GenseiGenerator(**self.config2)
        # self.gen = generators.gensei_wrapper.GenseiGenerator()
        logger.debug("1D structure generator started")
        # logger.debug("1D structure generator started")
        # import ipdb; ipdb.set_trace()
        self.gen.run()
        # logger.debug("vtk generated")
        # import ipdb; ipdb.set_trace()
        self._step1_deinit_save_stats(t0)
        logger.debug("step1 finished")


    def get_aposteriori_faces_and_vertices(self):
        """
        :return: (faces, vertices)
        """
        return self._aposteriori_surface_faces, self._aposteriori_surface_vertices

    def get_config_file_pattern(self):
        filepattern = self.config["filepattern"]
        filepattern = io3d.datawriter.filepattern_fill_slice_number_or_position(filepattern, "")
        base, ext = os.path.splitext(filepattern)
        return base + "_parameters.yaml"


    def __generate_vtk(self, vtk_file="~/tree.vtk"):
        logger.info("generating vtk for surface and volume compensated objects")
        vtk_file = op.expanduser(vtk_file)
        # from tree import TreeBuilder
        from .tb_vtk import TBVTK

        if "tree_data" in dir(self.gen):
            resolution = self.config["postprocessing"]["measurement_resolution"]
            method  = self.config["measurement"]["polygon_radius_selection_method"]
            tube_shape = self.config["measurement"]["tube_shape"]
            logger.debug("surface_tube_shape " + str(tube_shape))

            if method == "best":
                method_vol = "cylinder volume + sphere error"
                method_surf = "cylinder surface + sphere error + join error"
            else:
                method_vol = method
                method_surf = None

            # build volume tree
            logger.debug("vtk generation - volume compensated")
            tvg = TBVTK(
                cylinder_resolution=resolution,
                sphere_resolution=resolution,
                polygon_radius_selection_method=method_vol,
                tube_shape=tube_shape
            )
            # yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
            # tvg.importFromYaml(yaml_path)
            tvg.voxelsize_mm = self.voxelsize_mm
            tvg.shape = self.gen.areasize_px
            tvg.tube_skeleton = self.tube_skeleton
            output = tvg.buildTree()  # noqa
            # tvg.show()
            # TODO control output
            tvg.saveToFile(vtk_file)
            polydata_vol = tvg.polyData

            # build surface tree
            if method_surf is not None:
                logger.debug("vtk generation - surface compensated")
                from .tb_vtk import TBVTK
                tvg2 = TBVTK(
                    cylinder_resolution=resolution,
                    sphere_resolution=resolution,
                    polygon_radius_selection_method=method_surf,
                    tube_shape=tube_shape
                )
                tvg2.voxelsize_mm = self.voxelsize_mm
                tvg2.shape = self.gen.areasize_px
                tvg2.tube_skeleton = self.tube_skeleton
                output = tvg2.buildTree()  # noqa
                polydata_surf = tvg2.polyData
                # tvg.show()
                # TODO control output
                # tvg.saveToFile(vtk_file)
            else:
                polydata_surf = None

            return polydata_vol, polydata_surf

    def stop(self):
        self.stop_flag = True


    def filepattern_fill_potential_series(self):
        import io3d.datawriter
        # filepattern = self.config["filepattern"]
        filepattern = self.get_config_file_pattern()
        sn = io3d.datawriter.get_unoccupied_series_number(filepattern)
        filepattern = re.sub(r"({\s*})", r"", filepattern)
        filepattern = io3d.datawriter.filepattern_fill_series_number(filepattern, sn)
        return filepattern

    def filepattern_fill_series(self):
        """
        Return base and ext of file. The slice_number and slice_position is ignored.
        :return:
        """
        import io3d.datawriter
        filepattern = self.config["filepattern"]
        force_rewrite = self.config[CKEY_APPEARANCE]["force_rewrite"]
        # self.refresh_unoccupied_series_number()
        if force_rewrite and "filepattern_series_number" in self.config.keys():
            pass
        else:
            sn = io3d.datawriter.get_unoccupied_series_number(filepattern)
            self.config["filepattern_series_number"] = sn

        filepattern_series_number = self.config["filepattern_series_number"]

        filepattern = re.sub(r"({\s*})", r"", filepattern)

        filepattern = io3d.datawriter.filepattern_fill_series_number(filepattern, filepattern_series_number)
        return filepattern

    def filepattern_split(self):
        filepattern = self.filepattern_fill_series()
        filepattern = re.sub(r"({.*?})", r"", filepattern)
        root, ext = op.splitext(filepattern)
        return root, ext

    def get_fn_base(self):
        fn_base, fn_ext = self.filepattern_split()

        fn_base = op.expanduser(fn_base)
        return fn_base

    def refresh_unoccupied_series_number(self):
        config_filepattern = self.get_config_file_pattern()
        series_number = io3d.datawriter.get_unoccupied_series_number(filepattern=config_filepattern)
        # series_number = io3d.datawriter.get_unoccupied_series_number(filepattern=self.config["filepattern"])
        self.config['filepattern_series_number'] = series_number

    def save_parameters(self, filename=None):

        if filename is None:
            fn_base = self.get_fn_base()
            dirname = op.dirname(fn_base)
            if not op.exists(dirname):
                os.makedirs(dirname)
            filename = fn_base + "_config.yaml"

        io3d.misc.obj_to_file(self.config, filename=filename)
        return filename

    def save_config_and_measurement(self, filename=None):
        if filename is None:
            fn_base = self.get_fn_base()
            dirname = op.dirname(fn_base)
            if not op.exists(dirname):
                os.makedirs(dirname)
            filename = fn_base + "_config_and_measurement.yaml"
        logger.debug(f"save config and measurement to path: {filename}")
        config_and_measurement = self.get_config_and_measurement()

        io3d.misc.obj_to_file(config_and_measurement, filename=filename)

    def save_log(self):
        fn_base = self.get_fn_base()
        handler = logging.FileHandler(fn_base + ".log")
        handler.setFormatter(self.formatter)
        handler.setLevel(self.loglevel)
        self.memoryhandler.setTarget(handler)
        self.memoryhandler.flush()
        self.memoryhandler.flushLevel = self.loglevel

    def generate_volume(self):
        background_intensity=self.config["postprocessing"]["background_intensity"]
        self.data3d = self.gen.generate_volume(
            voxelsize_mm=self.config["areasampling"]["voxelsize_mm"],
            shape=self.config["areasampling"]["areasize_px"],
            tube_skeleton=self.tube_skeleton, dtype="uint8", background_intensity=background_intensity)
        # self.voxelsize_mm = self.gen.voxelsize_mm
        # this will change negative and positive
        id = self.config['generator_id']
        id = self.get_generator_id_by_name_or_number(id)
        # from pprint import pprint
        # plogger.debug(self.config)
        # import ipdb; ipdb.set_trace()
        # config = self._cfg_export_fcn[id](self.config)
        postprocessing_params = self.config["postprocessing"]
        postprocessing_params["negative"] = self._cfg_negative[id]
        data3d = self.postprocessing(**postprocessing_params)
        self.gen.data3d = data3d

    def step2(self):
        if self.parameters_changed_before_save:
            self.step1()
        # TODO split save_volume and save_parameters
        if len(self.tube_skeleton) == 0:
            logger.error("No data generated. 1D skeleton is empty.")
            return
        logger.debug(f"filepattern abspath: {self.config['filepattern_abspath']}")
        self.refresh_unoccupied_series_number()
        self.save_parameters()
        self.save_log()
        t0 = datetime.datetime.now()
        fn_base = self.get_fn_base()
        # config["filepattern"] = filepattern

        self._aposteriori_numeric_measurement(fn_base)
        logger.debug(f"step2 save stats: {fn_base}")
        self.save_stats(fn_base)
        t1 = datetime.datetime.now()
        self.save_1d_model_to_file(fn_base + "_vt.yaml")

        logger.debug("before volume generate " + str(t1 - t0))

        self.save_surface_to_file(fn_base + "_surface.vtk")
        #self.save_surface_to_file(fn_base + "_surface.vtk")
        # postprocessing
        skip_vg = self.config[CKEY_APPEARANCE]["skip_volume_generation"]

        t2 = t1
        if (not skip_vg) and ("generate_volume" in dir(self.gen)):
            # self.data3d = self.gen.generate_volume()
            self.generate_volume()
        # self.gen.saveVolumeToFile(self.config["filepattern"])

            t2 = datetime.datetime.now()
            logger.debug("volume generated in: " + str(t2 - t0))
            self.gen.saveVolumeToFile(self.config["filepattern_abspath"])
        t3 = datetime.datetime.now()
        logger.info("time before volume generate: " + str(t1 - t0))
        logger.info("time before volume save: " + str(t2 - t0))
        logger.info("time after volume save: " + str(t3 - t0))
        self.stats_times["step2_init_datetime"] = str(t0)
        self.stats_times["step2_numeric_measurement_time_s"] = (t1 - t0).total_seconds()
        self.stats_times["step2_generate_volume_time_s"] = (t2 - t1).total_seconds()
        self.stats_times["step2_save_volume_time_s"] = (t3 - t2).total_seconds()
        self.stats_times["step2_total_time_s"] = (t3 - t0).total_seconds()
        self.stats_times["step2_finish_datetime"] = str(t3)
        self.stats_times["step2_finished"] = True
        # fnp_abs = self.config["filepattern_abspath"]
        self.save_config_and_measurement(filename=fn_base + "_config_and_measurement.yaml")

        one_row_filename = self.config["output"]["one_row_filename"]
        logger.debug("write stats to common spreadsheet")
        if one_row_filename != "":
            # self.prepare_stats()
            self.save_stats_to_row(one_row_filename)
        else:
            self.prepare_stats()
        logger.debug("step2 finished")
        # self.memoryhandler.flush()

    def save_1d_model_to_file(self, outputfile):
        tree_data = dili.ndarray_to_list_in_structure(self.gen.tree_data)
        tree = {
            "voxelsize_mm": np.asarray(self.config["areasampling"]["voxelsize_mm"]).tolist(),
            "voxelsize_px": np.asarray(self.config["areasampling"]["areasize_px"]).tolist(),
            "Graph": {"0": tree_data}
        }
        io3d.misc.obj_to_file(tree, outputfile)


    def save_surface_to_file(self, outputfile, lc_all="C"):
        import vtk
        logger.debug("vtk version " + str(vtk.VTK_BUILD_VERSION))
        if lc_all is not None:
            import locale
            locale.setlocale(locale.LC_ALL, lc_all)
        # import ipdb; ipdb.set_trace()
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileName(outputfile)
        try:
            writer.SetInputData(self.polydata_volume)
        except:
            logger.warning("old vtk is used")
            writer.SetInput(self.polydata_volume)
        writer.Write()

    def postprocessing(
            self,
            gaussian_blur=True,
            gaussian_filter_sigma_mm=0.01,
            add_noise=False,
            # gaussian_noise_stddev=10.0,
            # gaussian_noise_center=0.0,
            limit_negative_intensities=True,
            noise_rng_seed=0,
            noise_exponent=0.0,
            noise_lambda0=0.02,
            noise_lambda1=1.0,
            noise_std=40.0,
            noise_mean=30.0,
            #            surface_measurement=False,
            #            measurement_multiplier=-1,
            measurement_resolution=20,
            output_dtype="uint8",
            intensity_profile_radius=[0.7, 1.0, 1.3],
            intensity_profile_intensity=[190, 200, 30],
            negative=False,
            background_intensity=20

    ):
        # negative is removed because we want it hide. The tab widget is used to control this
        # property
        dt = self.data3d.dtype

        logger.debug(f"len(unique(data3d)): {len(np.unique(self.data3d))}")
        logger.debug(f"describe(data3d) before gaussian: {scipy.stats.describe(self.data3d.flatten())}")
        if gaussian_blur:
            sigma_px = gaussian_filter_sigma_mm / self.voxelsize_mm

            self.data3d = scipy.ndimage.filters.gaussian_filter(
                self.data3d,
                sigma=sigma_px)
        logger.debug(f"describe(data3d) before noise: {scipy.stats.describe(self.data3d.flatten())}")

        if add_noise:
            noise = self.generate_noise()
            # noise = noise.astype(self.data3d.dtype)
            # noise = np.random.normal(loc=gaussian_noise_center, scale=gaussian_noise_stddev, size=self.data3d.shape)

            self.data3d = self.data3d + noise

        logger.debug(f"negative: {negative}")
        if negative:
            self.data3d = 255 - self.data3d
        logger.debug(f"after negative describe(data3d): {scipy.stats.describe(self.data3d.flatten())}")


        if limit_negative_intensities:
            #self.data3d[self.data3d < 0] = 0
            limit_ndarray(self.data3d, minimum=0, maximum=255)
        self.data3d = self.data3d.astype(dt)
        # self.config["postprocessing"]["measurement_multiplier"] = measurement_multiplier
        # negative = self.config["postprocessing"]["negative"] = measurement_multiplier

        return self.data3d

    def generate_noise(self):
        pparams = self.config["postprocessing"]
        parea = self.config["areasampling"]
        # data3d = self.postprocessing(**postprocessing_params)
        noise_params = dict(
            shape=parea["areasize_px"],
            sample_spacing=parea["voxelsize_mm"],
            exponent=pparams["noise_exponent"],
            random_generator_seed=pparams["noise_rng_seed"],
            lambda0=pparams["noise_lambda0"],
            lambda1=pparams["noise_lambda1"],
        )
        noise = ndnoise.noises(
            **noise_params
        ) #.astype(np.float16)
        noise = pparams["noise_std"] * noise / np.std(noise)
        noise = noise + pparams["noise_mean"]
        return noise

    def _aposteriori_numeric_measurement(self, fn_base):
        # import numpy as np
        from .tb_volume import TBVolume
        measurement_multiplier = self.config[CKEY_OUTPUT]["aposteriori_measurement_multiplier"]
        surface_measurement = self.config[CKEY_OUTPUT]["aposteriori_measurement"]

        vxsz = self.config["areasampling"]["voxelsize_mm"]
        vxsz = np.asarray(vxsz).astype(np.float) / measurement_multiplier
        shape = self.config["areasampling"]["areasize_px"]
        if measurement_multiplier > 0 and surface_measurement:
            shape = np.asarray(shape) * measurement_multiplier
            self._numeric_surface_measurement_shape = shape

            shape = shape.astype(np.int)
            tvgvol = TBVolume()
            tvgvol.voxelsize_mm = vxsz
            tvgvol.shape = shape
            tvgvol.tube_skeleton = self.gen.tree_data

            data3d = tvgvol.buildTree()
            import measurement
            surface, vertices, faces = measurement.surface_measurement(data3d, tvgvol.voxelsize_mm,
                                                                       return_vertices_and_faces=True)
            self._aposteriori_surface_vertices = vertices
            self._aposteriori_surface_faces = faces

            volume = np.sum(data3d > 0) * np.prod(vxsz)

            self.dataframes["overall"]["aposteriori numeric surface [mm^2]"] = [surface]
            self.dataframes["overall"]["aposteriori numeric volume [mm^3]"] = [volume]

            filename = fn_base + "_raw_{:06d}.jpg"
            import io3d.misc
            data = {
                'data3d': data3d.astype(np.uint8),  # * 70,  # * self.output_intensity,
                'voxelsize_mm': vxsz,
                # 'segmentation': np.zeros_like(self.data3d, dtype=np.int8)
            }
            io3d.write(data, filename)
        else:
            surface = 0
        return surface

    def get_stats(self):
        """ Return volume, surface, length and radius information.

        :return:
        Using one of generators to compute statistics.
        """
        from .generators import unconnected_cylinders as uncy
        gen = uncy.UnconnectedCylinderGenerator(
            areasize_px=self.config["areasampling"]["areasize_px"],
            voxelsize_mm=self.config["areasampling"]["voxelsize_mm"],
        )
        gen.init_stats()
        for id, element in self.tube_skeleton.items():
            if (
                                "nodeA_ZYX_mm" in element.keys() and
                                "nodeB_ZYX_mm" in element.keys() and
                                "radius_mm" in element.keys()
            ):
                nodeA = element["nodeA_ZYX_mm"]

                nodeB = element["nodeB_ZYX_mm"]
                radius = element["radius_mm"]
                gen.add_cylinder_to_stats(nodeA, nodeB, radius)
            else:
                logger.warning("Missing key in element id (" +
                               str(id) + "). Two point and radius are required")
        return gen.get_stats()

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

        # this compute correct even in case we are using cylinders
        # df = self.gen.get_stats()
        # this works fine even if data are loaded from file
        df = self.get_stats()
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
        dfoverallf["mean radius [mm]"] = df["radius"].mean()


        # surface and volume measurement
        import vtk
        mass = vtk.vtkMassProperties()
        # mass.SetInputData(object1Tri.GetOutput())
        mass.SetInputData(self.polydata_volume)
        vol = mass.GetVolume()
        if self.polydata_surface is None:
            surf = mass.GetSurfaceArea()
        else:
            mass = vtk.vtkMassProperties()
            mass.SetInputData(self.polydata_surface)
            surf = mass.GetSurfaceArea()
        dfoverallf["numeric volume [mm^3]"] = [vol]
        dfoverallf["numeric surface [mm^2]"] = [surf]
        dfoverallf["numeric volume fraction []"] = [vol/self.gen.area_volume]
        dfoverallf["negative numeric volume fraction []"] = [1. - vol/self.gen.area_volume]
        dfoverallf["negative numeric volume [mm^3]"] = [self.gen.area_volume - vol]
        self.dataframes["overall"] = dfoverallf

        st = self.stats_times
        # logger.debug("st ", st)
        note_df = pd.DataFrame([st], columns=st.keys())
        # logger.debug(note_df)
        # logger.debug(note_df.to_dict())

        self.dataframes["processing_info"] = note_df


    def load_tube_skeleton(self, filename):
        """ Load tube skeleton and remember it.

        :param filename:
        :return:
        """
        # params = io3d.misc.obj_from_file(filename=filename)
        from .import tree
        tube_skeleton, rawdata = tree.read_tube_skeleton_from_yaml(filename, return_rawdata=True)
        area = tree.parse_area_properties(rawdata)
        self.config["areasampling"].update(area)
        self.set_tube_skeleton(tube_skeleton)

    def set_tube_skeleton(self, tube_skeleton):
        self.tube_skeleton = tube_skeleton

    def get_tube_skeleton(self):
        return self.tube_skeleton

    def load_config(self, filename):
        """ Load config from file.
        :param filename:
        :return:
        """
        params = io3d.misc.obj_from_file(filename=filename)
        self.use_default_config()
        self.update_config(**params)

    def save_stats(self, fn_base):
        import pandas as pd

        for dfname in self.dataframes:
            df = self.dataframes[dfname]
            # to csv
            df.to_csv(fn_base + "_" + dfname + ".csv")

        try:
            writer = pd.ExcelWriter(fn_base + "_output.xlsx", engine="xlsxwriter")

            for dfname in ["overall", "density"]:
                logger.debug("adding xls list " + dfname)
                df = self.dataframes[dfname]
                # to excel
                df.to_excel(writer, dfname)
            writer.save()
        except:
            import traceback
            traceback.print_exc()
            s = traceback.format_exc()
            logger.warning(s)

    def get_flatten_config(self):
        """ Put input configuration into one row.

        :return:
        """
        config = self.config
        config_fl = dili.flatten_dict(config, join=lambda a, b: a + ' ' + b)
        config_fl = dict(config_fl)
        return config_fl

    def config_to_row_dataframe(self):
        config_fl = self.get_flatten_config()
        config_df = pd.DataFrame([config_fl], columns=config_fl.keys())
        return config_df

    def get_config_and_measurement(self):
        self.prepare_stats()
        dfo = self.dataframes["overall"].to_dict(orient="records")[0]
        dfd = self.dataframes["density"].to_dict(orient="records")[0]
        dfi = self.dataframes["processing_info"].to_dict(orient="records")[0]

        dfo.update(dfd)
        data_structured = {
            "measurement": dfo,
            "processing_info": dfi,
            "config": self.config
        }
        return data_structured

    def get_flatten_config_and_measurement(self):
        import imtools.dili
        return imtools.dili.flatten_dict_join_keys(self.get_config_and_measurement(), " ")

    def save_stats_to_row(self, filename, note=""):
        """ Save stats to row

        :param filename:
        :param note:
        :return:
        """
        self.prepare_stats()
        filename = op.expanduser(filename)
        import pandas as pd
        # filename = op.expanduser(filename)
        # dfo = self.dataframes["overall"]
        # dfd = self.dataframes["density"]
        # dfi = self.dataframes["processing_info"]
        #
        #
        # # values must be a list for dataframe
        # # new_values = []
        # # for val in config_fl.values():
        # #     new_values.append([val])
        #
        # # config_fl_li = dict(zip(config_fl.keys(), new_values))
        # # config_df = pd.DataFrame(config_fl_li)
        # config_fl = self.get_flatten_config()
        #
        # config_df = pd.DataFrame([config_fl], columns=config_fl.keys())
        # # import ipdb; ipdb.set_trace()
        # dfout = pd.concat([dfi, dfo, dfd, config_df], axis=1)
        config_fl = self.get_flatten_config_and_measurement()

        dfout = pd.DataFrame([config_fl], columns=config_fl.keys())
        if op.exists(filename):
            dfin = pd.read_csv(filename)
            dfout = pd.concat([dfin, dfout], axis=0, sort=False)
        else:
            dirname = op.dirname(filename)
            if not op.exists(dirname):
                os.makedirs(dirname)

        dfout.to_csv(filename, index=False)
        # import ipdb; ipdb.set_trace()
        # pass

    def _config2generator_general_export(self, config):
        return {
            'voxelsize_mm': config["areasampling"]["voxelsize_mm"],
            'areasize_px': config["areasampling"]["areasize_px"],
            "intensity_profile_radius": config["postprocessing"]["intensity_profile_radius"],
            "intensity_profile_intensity": config["postprocessing"]["intensity_profile_intensity"],
            "tube_shape": config["measurement"]["tube_shape"]
        }

    def _config2generator_tubes_export(self, config):
        config["postprocessing"]["negative"] = False
        generator_config = self._config2generator_general_export(config)
        return generator_config

    def _config2generator_porosity_export(self, config):
        config["postprocessing"]["negative"] = True
        generator_config = self._config2generator_general_export(config)
        return generator_config

    def _config2generator_gensei_export(self, config):
        asp = config["areasampling"]
        vs_mm = np.asarray(asp["voxelsize_mm"])
        resolution = 1.0 / vs_mm
        dct = {
            'dims': asp["areasize_mm"],
            'n_slice': asp["areasize_px"][0],
            'resolution': [resolution[1], resolution[2]]
        }
        return dct

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
            self.step1()
            self.step2()


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
        self.logfile = op.join(self.config_dir, log_file)

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


# @click.command()
# @begin.start
def new_main(
        parameterfile=None,
        debug=True,
        d=True,
        nointeractivity=False,
        logfile="~/teigen.log",
):
    """ Run test image generator.

    :param parameterfile:
    :param debug:
    :param d:
    :param nointeractivity:
    :param logfile:
    :return:
    """
    logger = logging.getLogger()

    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    logger.addHandler(ch)

    config_file_manager = ConfigFileManager("teigen")
    config_file_manager.init_config_dir()

    if parameterfile is None:
        parameterfile = config_file_manager.init_config_file

    if debug or d:
        ch.setLevel(logging.DEBUG)

    # default param file
    if not op.exists(op.expanduser(parameterfile)):
        parameterfile = None

    if nointeractivity:
        tg = Teigen(logfile=logfile)
        if parameterfile is not None:
            params = io3d.misc.obj_from_file(parameterfile)
            tg.update_config(**params)
        tg.step1()
        # tg.run(**params)
        tg.step2()
    else:
        from PyQt5.QtWidgets import QApplication
        from .gui import TeigenWidget
        app = QApplication(sys.argv)
        params = None
        if parameterfile is not None:
            params = io3d.misc.obj_from_file(parameterfile)
        cw = TeigenWidget(logfile=logfile, config=params)
        cw.show()
        app.exec_()

def limit_ndarray(data, minimum, maximum):
    data[data < minimum] = minimum
    data[data > maximum] = maximum
    return data

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

    # default param file
    if not op.exists(op.expanduser(args.parameterfile)):
        args.parameterfile = None

    if args.nointeractivity:
        tg = Teigen(logfile=args.logfile)
        if args.parameterfile is not None:
            params = io3d.misc.obj_from_file(args.parameterfile)
            tg.update_config(**params)
        tg.step1()
        # tg.run(**params)
        tg.step2()
    else:
        from PyQt5.QtWidgets import QApplication
        from .gui import TeigenWidget
        app = QApplication(sys.argv)
        params = None
        if args.parameterfile is not None:
            try:
                params = io3d.misc.obj_from_file(args.parameterfile)
            except Exception as e:
                import traceback
                logger.warning(f"Problem with reading: {args.parameterfile}")
                logger.warning(traceback.format_exc())
        cw = TeigenWidget(logfile=args.logfile, config=params)
        cw.show()
        app.exec_()
