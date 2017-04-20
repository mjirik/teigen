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

import generators.cylinders
import generators.gensei_wrapper
import generators.unconnected_cylinders
import dili
import io3d.datawriter
import io3d.misc
import dictwidgetqt

CKEY_APPEARANCE = "appearance"
CKEY_OUTPUT = "output"

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
        self.version = "0.2.17"
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
        self.fig_3d_render_snapshot = None

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
        config["filepattern"] = "~/teigen_data/{seriesn:03d}/data{:06d}.jpg"
        # config['filepattern_series_number'] = series_number
        # self.config["voxelsize_mm"] = [1.0, 1.0, 1.0]
        # self.config["areasize_mm"] = [100.0, 100.0, 100.0]
        # self.config["areasize_px"] = [100, 100, 100]
        config[CKEY_APPEARANCE] = {
            "show_aposteriori_surface": True,
            "skip_volume_generation": False,
            "noise_preview": False,
        }
        config["output"] = {
            "one_row_filename": "~/teigen_data/output_rows.csv",
            "aposteriori_measurement": False,
            "aposteriori_measurement_multiplier": 1.0,
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
        logger.debug("1D structure generator started")
        #print "1D structure generator started"
        #import ipdb; ipdb.set_trace()
        self.gen.run()

        logger.debug("1D structure is generated")
        self.polydata = self.__generate_vtk(self.temp_vtk_file)
        #logger.debug("vtk generated")
        #import ipdb; ipdb.set_trace()

        self.prepare_stats()
        one_row_filename = self.config["output"]["one_row_filename"]
        if one_row_filename != "":
            self.save_stats_to_row(one_row_filename)



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
        # self.refresh_unoccupied_series_number()
        if "filepattern_series_number" not in self.config.keys():
            sn = io3d.datawriter.get_unoccupied_series_number(filepattern)
            self.config["filepattern_series_number"] = sn

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
        # if self.parameters_changed_before_save:
        #     self.run()
        # prepare path to save


        fn_base = None
        if filename is None:
            self.refresh_unoccupied_series_number()
            # config_filepattern = self.get_config_file_pattern()
            # series_number = io3d.datawriter.get_unoccupied_series_number(filepattern=config_filepattern)
            # # series_number = io3d.datawriter.get_unoccupied_series_number(filepattern=self.config["filepattern"])
            # self.config['filepattern_series_number'] = series_number

            fn_base = self.get_fn_base()

            dirname = op.dirname(fn_base)
            if not op.exists(dirname):
                os.makedirs(dirname)
            filename = fn_base + "_parameters.yaml"


        io3d.misc.obj_to_file(self.config, filename=filename)

    def save_log(self):
        fn_base = self.get_fn_base()
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
        self.save_log()
        import io3d.misc
        import time
        t0 = time.time()
        fn_base = self.get_fn_base()
        # config["filepattern"] = filepattern

        self._numeric_measurement(fn_base)
        self.save_stats(fn_base)
        t1 = time.time()
        logger.debug("before volume generate " + str(t1-t0))
        # postprocessing
        skip_vg = self.config[CKEY_APPEARANCE]["skip_volume_generation"]
        if (not skip_vg) and ("generate_volume" in dir(self.gen)):
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
            # gaussian_noise_stddev=10.0,
            # gaussian_noise_center=0.0,
            limit_negative_intensities=True,
            noise_rng_seed=0,
            noise_exponent=0.0001,
            noise_lambda_start=0.1,
            noise_lambda_stop=3.0,
            noise_amplitude = 40.0,
            noise_mean = 30.0,
            #            surface_measurement=False,
            #            measurement_multiplier=-1,
            measurement_resolution=20,
            output_dtype="uint8",
            negative=False,


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

        if negative:
            self.data3d = 255 - self.data3d

        if limit_negative_intensities:
            self.data3d[self.data3d < 0] = 0
        # self.config["postprocessing"]["measurement_multiplier"] = measurement_multiplier
        # negative = self.config["postprocessing"]["negative"] = measurement_multiplier

        return self.data3d

    def generate_noise(self):
        import ndnoise
        import ndnoise.generator
        pparams = self.config["postprocessing"]
        # data3d = self.postprocessing(**postprocessing_params)
        noise = ndnoise.noises(
            shape=self.gen.areasize_px,
            sample_spacing=self.gen.voxelsize_mm,
            exponent=pparams["noise_exponent"],
            random_generator_seed=pparams["noise_rng_seed"],
            lambda_start=pparams["noise_lambda_start"],
            lambda_stop=pparams["noise_lambda_stop"],

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
        measurement_multiplier = self.config[CKEY_OUTPUT]["aposteriori_measurement_multiplier"]
        surface_measurement = self.config[CKEY_OUTPUT]["aposteriori_measurement"]

        vxsz = self.config["areasampling"]["voxelsize_mm"]
        vxsz = np.asarray(vxsz).astype(np.float) / measurement_multiplier
        shape = self.config["areasampling"]["areasize_px"]
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
            self._aposteriori_surface_vertices = vertices
            self._aposteriori_surface_faces = faces

            volume = np.sum(data3d > 0) * np.prod(vxsz)

            self.dataframes["overall"]["aposteriori numeric surface [mm^2]"] = [surface]
            self.dataframes["overall"]["aposteriori numeric volume [mm^3]"] = [volume]

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


        # surface and volume measurement
        import vtk
        mass = vtk.vtkMassProperties()
        # mass.SetInputData(object1Tri.GetOutput())
        mass.SetInputData(self.polydata)
        surf = mass.GetSurfaceArea()
        vol = mass.GetVolume()
        dfoverallf["numeric volume [mm^3]"]=[vol]
        dfoverallf["numeric surface [mm^2]"]=[surf]
        self.dataframes["overall"] = dfoverallf

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

    def save_stats_to_row(self, filename, note=""):
        import pandas as pd
        filename = op.expanduser(filename)
        tm = datetime.datetime.now().isoformat()
        note_df = pd.DataFrame({"datetime": [tm], "note": [note]})
        dfo = self.dataframes["overall"]
        dfd = self.dataframes["density"]



        config = self.config
        config_fl = dili.flattenDict(config, join=lambda a,b:a+' '+b)
        config_fl = dict(config_fl)

        # values must be a list for dataframe
        new_values = []
        for val in config_fl.values():
            new_values.append([val])

        config_fl_li = dict(zip(config_fl.keys(), new_values))
        config_df = pd.DataFrame(config_fl_li)
        # import ipdb; ipdb.set_trace()
        dfout = pd.concat([note_df, dfo, dfd, config_df], axis=1)

        if op.exists(filename):

            dfin = pd.read_csv(filename)
            dfout = pd.concat([dfin, dfout], axis=0)

        dfout.to_csv(filename, index=False)
        # import ipdb; ipdb.set_trace()
        # pass


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
        from PyQt4.QtGui import QApplication
        from gui import TeigenWidget
        app = QApplication(sys.argv)
        params = None
        if args.parameterfile is not None:
            params = io3d.misc.obj_from_file(args.parameterfile)
        cw = TeigenWidget(logfile=args.logfile, config=params)
        cw.show()
        app.exec_()
