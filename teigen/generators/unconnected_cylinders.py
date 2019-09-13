#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © %YEAR% %USER% <%MAIL%>
#
# Distributed under terms of the %LICENSE% license.

"""
%HERE%
"""

import logging
logger = logging.getLogger(__name__)
import argparse
import numpy as np
# import ..geometry3 as g3
# from ..geometry3d import plane_fit
from .. import geometry3d as g3
import os.path
from . import general


def __half_plane(self, perp, plane_point, point):
    cdf = (np.array(point) - np.array(plane_point))
    out = perp[0] * cdf[0] + \
          perp[1] * cdf[1] + \
          perp[2] * cdf[2]
    return out > 0


def _const(value):
    return value


class UnconnectedCylinderGenerator(general.GeneralGenerator):
    def __init__(self,
                 build=True,
                 gtree=None,
                 # endDistMultiplicator=1,
                 # use_joints=True,
                 voxelsize_mm=[1.0, 1.0, 1.0],
                 # voxelsize_mm_z=1.0,
                 # voxelsize_mm_x=1.0,
                 # voxelsize_mm_y=1.0,
                 areasize_px=[100, 100, 100],
                 # area_shape_z=100,
                 # area_shape_x=100,
                 # area_shape_y=100,
                 element_number=-1,
                 radius_distribution="normal",
                 radius_distribution_uniform=None,
                 radius_distribution_normal=None,
                 radius_distribution_fixed=None,
                 radius_distribution_minimum=0.05,
                 radius_distribution_maximum=10.0,
                 radius_distribution_mean=0.1,
                 radius_distribution_standard_deviation=0.1,
                 length_distribution_mean=1.0,
                 length_distribution_standard_deviation=0.1,
                 # intensity_profile=None
                 intensity_profile_radius=[0.4, 0.7, 1.0, 1.3],
                 intensity_profile_intensity=[195, 190, 200, 30],
                 orientation_anisotropic=True,
                 orientation_alpha_rad=0.0,
                 orientation_beta_rad=0.0,
                 orientation_variance_rad=0.1,
                 allow_overlap=False,
                 volume_fraction=0.1,
                 maximum_1000_iteration_number=10,
                 random_generator_seed=0,
                 last_element_can_be_smaller=False,
                 tube_shape=True
                 ):
        """
        gtree is information about input data structure.
        endDistMultiplicator: make cylinder shorter by multiplication of radius
        intensity_profile: Dictionary type. Key is radius and value is required intensity.
        @param tube_shape: create tube shape if true, otherwise create cylinders
        """
        # area_shape = [area_shape_z,area_shape_x, area_shape_y]
        # voxelsize_mm = [
        #     voxelsize_mm_z,
        #     voxelsize_mm_x,
        #     voxelsize_mm_y
        # ]
        super(general.GeneralGenerator, self).__init__()
        # general.GeneralGenerator.__init__(self)
        self.build = build
        # self.filename = "output{:05d}.jpg"
        self.areasize_px = np.asarray(areasize_px)
        self.voxelsize_mm = np.asarray(voxelsize_mm)
        self.element_number = element_number
        self.radius_maximum = radius_distribution_maximum
        self.radius_minimum = radius_distribution_minimum
        # self.intensity_profile = intensity_profile
        self.intensity_profile = dict(zip(intensity_profile_radius, intensity_profile_intensity))
        self._cylinder_nodes = []
        self._cylinder_nodes_radiuses = []
        self.random_generator_seed = random_generator_seed
        self.radius_generator = _const
        self.radius_generator_args = [radius_distribution_mean]
        self.area_volume = np.prod(self.areasize_px * self.voxelsize_mm)
        if radius_distribution_normal is not None:
            logger.warning("Deprecated use of radius_distribution_normal. Use radius_distribution='normal'")
        else:
            if radius_distribution == "normal":
                radius_distribution_normal = True
                radius_distribution_fixed = False
                radius_distribution_uniform = False
            elif radius_distribution == "fixed":
                radius_distribution_normal = False
                radius_distribution_fixed = True
                radius_distribution_uniform = False
            elif radius_distribution == "uniform":
                radius_distribution_normal = False
                radius_distribution_fixed= False
                radius_distribution_uniform = True
        if radius_distribution_uniform:
            self.radius_generator = np.random.uniform
            self.radius_generator_args = [radius_distribution_minimum, radius_distribution_maximum]
        if radius_distribution_normal:
            self.radius_generator = general.random_normal
            self.radius_generator_args = [radius_distribution_mean, radius_distribution_standard_deviation]
        self.alow_overlap = allow_overlap

        self.length_generator = general.random_normal
        self.length_generator_args = [length_distribution_mean, length_distribution_standard_deviation]

        self.requeseted_volume_fraction = volume_fraction
        self.max_iteration = 1000 * maximum_1000_iteration_number
        self.last_element_can_be_smaller = last_element_can_be_smaller
        # import ipdb; ipdb.set_trace()
        # input of geometry and topology
        # self.V = []
        # self.CV = []
        # self.joints = {}
        # self.joints_lar = []
        # self.gtree = gtree
        # self.endDistMultiplicator = endDistMultiplicator
        # self.use_joints = use_joints
        self.surface = 0
        # self.LEN_STEP_CONSTANT = 0.1
        self.MAKE_IT_SHORTER_CONSTANT = 2.0
        # self.DIST_MAX_RADIUS_MULTIPLICATOR = 3.0
        self.OVERLAPS_ALOWED = False
        self.tree_data = {}
        self.progress_callback = None
        # self.collision_model = g3.CollisionModelSpheres(areasize=(self.areasize_px * self.voxelsize_mm))
        self.collision_model = g3.CollisionModelCombined(areasize=(self.areasize_px * self.voxelsize_mm))

        self.area_volume = np.prod(self.voxelsize_mm * self.areasize_px)
        self.orientation_anisotropic = orientation_anisotropic
        self.orientation_alpha_rad = orientation_alpha_rad
        self.orientation_beta_rad = orientation_beta_rad
        self.orientation_variance_rad = orientation_variance_rad
        self.tube_shape = tube_shape

    def _add_cylinder_if_no_collision(self, pt1, pt2, radius,
                                      COLLISION_RADIUS=1.5  # higher then sqrt(2)
                                      ):
        if self.alow_overlap:
            return self.collision_model.add_tube(pt1, pt2, radius)
        else:
            return self.collision_model.add_tube_if_no_collision(pt1, pt2, radius)

    def run(self):
        logger.info("cylynder generator running")

        self.tree_data = {

        }
        np.random.seed(self.random_generator_seed)
        self.surface = 0
        # pts = np.random.random([self.element_number, 3]) * self.areasize_px * self.voxelsize_mm

        # construct voronoi
        import scipy.spatial
        import itertools
        self.init_stats()

        # radius = self.radius_maximum
        # for i, two_points in enumerate(vor3.ridge_points):
        # for i in range(self.element_number):
        self.generation_break_causes = {
            "radius_maximum": 0,
            "radius_minimum": 0,
            "collision": 0,
            "radius_bigger_than_areasize": 0,
            "length_bigger_than_areasize": 0,


        }
        while not self.is_final_iteration():
            self.create_cylinder()
        logger.info(self.generation_break_causes)
        self.get_stats()
        logger.debug(self.generation_break_causes)
        self.data3d = None

    def is_final_iteration(self):
        if self.stop:
            return True
        self.iterations += 1
        stats = self.get_stats()

        if self.element_number > 0:
            n = len(self.geometry_data["volume"])
            if n >= self.element_number:
                return True
        self.actual_object_volume = np.sum(self.geometry_data["volume"])
        actual_volume_fraction = self.actual_object_volume / self.area_volume

        logger.debug("iteration: " + str(self.iterations) + " / " + str(self.max_iteration))
        logger.debug("actual_volume_fraction: " + str(actual_volume_fraction))
        if self.iterations > self.max_iteration:
            return True
        elif actual_volume_fraction > self.requeseted_volume_fraction:
            return True
        else:
            return False

    def init_stats(self):
        self.stop = False
        self.iterations = 0
        self.actual_object_volume = 0
        self.requeseted_volume = self.requeseted_volume_fraction * self.area_volume
        self.geometry_data = {
            "length": [],
            "radius": [],
            "volume": [],
            "surface": [],
            "vector": [],
            "point1": [],
            "point2": []
        }

    def add_cylinder_to_stats(self, pt1, pt2, radius):
        pt1 = np.asarray(pt1)
        pt2 = np.asarray(pt2)
        edge = {
            "nodeA_ZYX_mm": pt1,
            "nodeB_ZYX_mm": pt2,
            # "radius_mm": radius
            # "radius_mm": 1 + np.random.rand() * (self.max_radius -1 )
            "radius_mm": radius,
        }
        self.tree_data[len(self.tree_data)] = edge
        # line_nodes = g3.get_points_in_line_segment(pt1, pt2, radius)
        # self._cylinder_nodes.extend(line_nodes)
        length = np.linalg.norm(pt1 - pt2)

        # if it is tube (pill)
        if self.tube_shape:
            surf = g3.tube_surface(radius, length)
            volume = g3.tube_volume(radius, length)
        else:
            # it is cylinder
            surf = 2 * np.pi * radius * (radius + length)
            volume =  np.pi * radius**2 * length

        vector = pt1 - pt2

        # TODO rename - add units
        self.geometry_data["length"].append(length)
        self.geometry_data["surface"].append(surf)
        self.geometry_data["radius"].append(radius)
        self.geometry_data["volume"].append(volume)
        self.geometry_data["vector"].append(vector)
        self.geometry_data["point1"].append(pt1)
        self.geometry_data["point2"].append(pt2)
        # self.geometry_data["collide_with_prevs"].append(collide_with_prevs)
        self.surface += surf

    def create_cylinder(
            self,
            try_shorter_iteration_number=8,
            n_nearest=4,
            length_to_radius_ratio=4
    ):
        generated = False
        while not generated:
            self.iterations += 1
            if self.is_final_iteration():
                if self.progress_callback is not None:
                    self.progress_callback(self, 1.0, statusbar_text="Skeleton created")
                return
            progress = self.iterations / (1. * self.max_iteration)
            # logger.debug("progress " + str(progress))
            object_volume = np.sum(self.geometry_data["volume"])
            actual_volume_fraction = object_volume / self.area_volume

            statusbar_text = str(self.iterations) + " Vv " + str(actual_volume_fraction)
            if self.progress_callback is not None:
                self.progress_callback(self, progress, statusbar_text=statusbar_text)
            # logger.debug(progress)
            # logger.debug(self.radius_generator_args)
            radius = self.radius_generator(*self.radius_generator_args)
            if radius > self.radius_maximum:
                self.generation_break_causes["radius_maximum"] += 1
                continue
            if radius < self.radius_minimum:
                self.generation_break_causes["radius_minimum"] += 1
                continue
            if (radius > (self.areasize_px * self.voxelsize_mm)).all():
                self.generation_break_causes["radius_bigger_than_areasize"] += 1
                continue
            # pt1 = self.collision_model.get_random_point(radius=radius)
            # pt2 = self.collision_model.get_random_point(radius=radius)
            # pt1 = np.random.random([3]) * self.areasize_px * self.voxelsize_mm
            # pt2 = np.random.random([3]) * self.areasize_px * self.voxelsize_mm
            # pt1, pt2 = self._make_cylinder_shorter(pt1, pt2, radius*self.MAKE_IT_SHORTER_CONSTANT)


            center = np.random.random([3]) * self.areasize_px * self.voxelsize_mm
            if self.collision_model.object_number > 2 * n_nearest:
                if self.iterations % 5:
                    # if self.collision_model.get_node_number() > n_nearest:
                    center = np.asarray(center)
                    npts, indexes, lengths = self.collision_model.n_closest_end_points(center, n_nearest)
                    center = np.mean(npts, axis=0)

            if self.orientation_anisotropic:
                # direction_vector = np.asarray(self.orientation_main).reshape(3, 1)
                # past solution
                # direction_vector = np.random.normal(direction_vector, self.orientation_variance_rad)
                # direction_vector = direction_vector / np.linalg.norm(direction_vector)
                direction_vector = g3.random_vector_along_direction(
                    alpha=self.orientation_alpha_rad, beta=self.orientation_beta_rad,
                    sigma=self.orientation_variance_rad, size=1
                )
                direction_vector = direction_vector.squeeze()

            else:
                direction_vector = g3.random_direction_vector()
            # direction_vector = np.asarray([0, 2**-0.5, 2**-0.5])
            # direction_vector = np.asarray([0, 2, 0])
            length = self.length_generator(*self.length_generator_args)
            if (length > (self.areasize_px * self.voxelsize_mm)).all():
                self.generation_break_causes["length_bigger_than_areasize"] += 1
                continue

            if self.tube_shape:
                volume = g3.tube_volume(radius, length)
            else:
                volume = g3.cylinder_volume(radius, length)

            planned_volume_is_too_much = \
                (( self.actual_object_volume + volume) / self.area_volume) > self.requeseted_volume_fraction
            if planned_volume_is_too_much and self.last_element_can_be_smaller:
                # just in case of last element and if is this feature enabled
                radius, length = self.pill_parameter_suggestion_for_last_object(radius, length)
                # pokud je navrhovaný objem přílišný
                pt1 = np.asarray(g3.translate(center, direction_vector, 0.5 * length))
                pt2 = np.asarray(g3.translate(center, direction_vector, -0.5 * length))
                collision = self._add_cylinder_if_no_collision(pt1, pt2, radius)

            else:
                # normal run

                pt1 = np.asarray(g3.translate(center, direction_vector, 0.5 * length))
                pt2 = np.asarray(g3.translate(center, direction_vector, -0.5 * length))

                try_shorter_i = 0
                collision = self._add_cylinder_if_no_collision(pt1, pt2, radius)
                while (collision is True and try_shorter_i < try_shorter_iteration_number):
                    try_shorter_i += 1
                    pt1, pt2 = g3.get_points_closer(pt1, pt2, relative_length=0.75)
                    collision = self._add_cylinder_if_no_collision(pt1, pt2, radius)

            if not collision:
                generated = True
            else:
                self.generation_break_causes["collision"] += 1

        if generated:
            self.add_cylinder_to_stats(pt1, pt2, radius=radius)
        # else:
        #     logger.debug(self.generation_break_causes)
        return

    def pill_parameter_suggestion_for_last_object(self, first_radius, first_length):
        length = first_length
        radius = g3.tube_radius_from_volume(self.actual_object_volume - self.requeseted_volume, length)

        # for alpha in np.linspace(1., 0., 10):
        #     length = alpha * first_length
        #     radius = g3.pill_radius_from_volume(self.actual_object_volume - self.requeseted_volume, length)
        #
        #     if radius >= self.radius_minimum:
        #         break


        if radius < self.radius_minimum:
            radius = self.radius_minimum

        return radius, length

    def get_stats(self):
        # self.assertTrue(False)
        import pandas as pd
        df = pd.DataFrame(self.geometry_data)
        # desc = df.describe()

        return df

    def _make_cylinder_shorter(self, nodeA, nodeB, radius):  # , radius, cylinder_id):
        vector = (np.asarray(nodeA) - np.asarray(nodeB)).tolist()
        if np.linalg.norm(vector) < 2 * radius:
            return None, None

        # mov circles to center of cylinder by size of radius because of joint
        nodeA = g3.translate(nodeA, vector,
                             -radius)  # * self.endDistMultiplicator)
        nodeB = g3.translate(nodeB, vector,
                             radius)  # * self.endDistMultiplicator)
        return nodeA, nodeB

    def _is_in_area(self, node, radius=None):
        return self.collision_model.is_point_in_area(node, radius)

    def add_cylinder(self, nodeA, nodeB, radius, cylinder_id):

        try:
            idA = tuple(nodeA)  # self.gtree.tree_data[cylinder_id]['nodeIdA']
            idB = tuple(nodeB)  # self.gtree.tree_data[cylinder_id]['nodeIdB']
        except:
            idA = 0
            idB = 0
            self.use_joints = False

        # vect = nodeA - nodeB
        # self.__draw_circle(nodeB, vect, radius)

        vector = (np.asarray(nodeA) - np.asarray(nodeB)).tolist()

        # mov circles to center of cylinder by size of radius because of joint
        nodeA = g3.translate(nodeA, vector,
                             -radius * self.endDistMultiplicator)
        nodeB = g3.translate(nodeB, vector,
                             radius * self.endDistMultiplicator)

        if all(nodeA == nodeB):
            logger.error("End points are on same place")

        ptsA, ptsB = g3.cylinder_circles(nodeA, nodeB, radius,
                                         element_number=30)
        CVlistA = self.__construct_cylinder_end(ptsA, idA)
        CVlistB = self.__construct_cylinder_end(ptsB, idB)

        CVlist = CVlistA + CVlistB

        self.CV.append(CVlist)


# lar add ball
#         ball0 = mapper.larBall(radius, angle1=PI, angle2=2*PI)([10, 16])
#         V, CV = ball0
#         # mapper.T
#         # ball = STRUCT(MKPOLS(ball0))
#
#         # mapper.T(1)(nodeA[0])(mapper.T(2)(nodeA[1])(mapper.T(3)(nodeA[1])(ball)))
#
#         lenV = len(self.V)
#
#         self.V = self.V + (np.array(V) + np.array(nodeA)).tolist()
#         self.CV = self.CV + (np.array(CV) + lenV).tolist()


def main():
    logging.basicConfig()
    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)

    # input parser
    parser = argparse.ArgumentParser(
        description='Histology analyser reporter. Try: \
python src/tb_volume.py -i ./tests/hist_stats_test.yaml'
    )
    parser.add_argument(
        '-i', '--inputfile',
        default=None,
        # required=True,
        help='input file, yaml file'
    )
    parser.add_argument(
        '-o', '--outputfile',
        default=None,
        help='output file, .raw, .dcm, .tiff, given by extension '
    )
    parser.add_argument(
        '-ot', '--outputfiletype',
        default='pkl',
        help='output file type.  raw, dcm, tiff, or pkl,   default is pkl, '
    )
    parser.add_argument(
        '-vs', '--voxelsize',
        default=[1.0, 1.0, 1.0],
        type=float,
        metavar='N',
        nargs='+',
        help='size of voxel (ZYX)'
    )
    parser.add_argument(
        '-ds', '--datashape',
        default=[200, 200, 200],
        type=int,
        metavar='N',
        nargs='+',
        help='size of output data in pixels for each axis (ZYX)'
    )
    parser.add_argument(
        '-g', '--generator',
        default='vol',
        type=str,
        help='Volume or surface model can be generated by use this option. \
                Use "vol", "volume" for volumetric model. For LAR surface model\
                use "lar". For VTK file use "vtk".'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='Debug mode')
    parser.add_argument(
        '-l', '--useLar', action='store_true',
        help='Use LAR')
    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)

    # startTime = datetime.now()

    generator_params = None
    generator_class = args.generator

    # if args.generator == "vtk":
    #     import gen_vtk_tree
    #     gen_vtk_tree.vt2vtk_file(args.inputfile, args.outputfile)
    #     return
    cylgen = CylinderGenerator()
    cylgen.step1()

    # tg = TreeBuilder(generator_class, generator_params)
    # tg.importFromYaml(args.inputfile)
    # tg.voxelsize_mm = args.voxelsize
    # tg.shape = args.datashape
    # tg.use_lar = args.useLar
    # data3d = tg.buildTree()
    #
    # logger.info("TimeUsed:" + str(datetime.now() - startTime))
    # # volume_px = sum(sum(sum(data3d)))
    # # volume_mm3 = volume_px * \
    # #     (tg.voxelsize_mm[0] * tg.voxelsize_mm[1] * tg.voxelsize_mm[2])
    # # logger.info("Volume px:" + str(volume_px))
    # # logger.info("Volume mm3:" + str(volume_mm3))
    #
    # # vizualizace
    # logger.debug("before visualization")
    # tg.show()
    # logger.debug("after visualization")

    # ukládání do souboru
    # if args.outputfile is not None:
    #     tg.saveToFile(args.outputfile, args.outputfiletype)


# class TreeGenerator(TreeConstructor):
#     """
#     back compatibility
#     """
#     pass

if __name__ == "__main__":
    main()
