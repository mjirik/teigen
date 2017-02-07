# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)
import unittest
import numpy as np
import teigen.geometry3d as g3

class GeometryTestCase(unittest.TestCase):
    def test_collision_model_no_collision(self):
        cm = g3.CollisionModelSpheres(areasize=[150, 151, 155])
        distance = 30
        pt1 = [20, 20, 20]
        pt2 = [20, 20, 60]
        pt3 = [20, 20 + distance, 20]
        pt4 = [20, 20 + distance, 60]
        r1 = 10
        r2 = 10
        collision1 = cm.add_cylinder_if_no_collision(pt1, pt2, radius=r1)
        collision2 = cm.add_cylinder_if_no_collision(pt3, pt4, radius=r2)

        self.assertEqual(collision2, False)


    def test_collision_model_collision(self):
        cm = g3.CollisionModelSpheres(areasize=[150, 151, 155])
        distance = 25
        pt1 = [20, 20, 20]
        pt2 = [20, 20, 60]
        pt3 = [20, 20 + distance, 20]
        pt4 = [20, 20 + distance, 60]
        r1 = 10
        r2 = 10
        collision1 = cm.add_cylinder_if_no_collision(pt1, pt2, radius=r1)
        collision2 = cm.add_cylinder_if_no_collision(pt3, pt4, radius=r2)

        self.assertEqual(collision2, True)

    def test_collision_model_out_of_area(self):
        """
        Cylinder end point should be in safe distance from boundary
        :return:
        """
        cm = g3.CollisionModelSpheres(areasize=[150, 151, 65])
        distance = 25
        pt1 = [20, 20, 20]
        pt2 = [20, 20, 60]
        pt3 = [20, 20 + distance, 20]
        pt4 = [20, 20 + distance, 60]
        r1 = 10
        r2 = 10
        collision1 = cm.add_cylinder_if_no_collision(pt1, pt2, radius=r1)
        collision2 = cm.add_cylinder_if_no_collision(pt3, pt4, radius=r2)

        self.assertEqual(collision2, True)

    def test_collision_model_n_closest(self):
        cm = g3.CollisionModelSpheres(areasize=[150, 151, 155])
        distance = 30
        pt1 = [20, 20, 20]
        pt2 = [20, 20, 60]
        pt3 = [20, 20 + distance, 20]
        pt4 = [20, 20 + distance, 60]
        r1 = 10
        r2 = 10
        collision1 = cm.add_cylinder_if_no_collision(pt1, pt2, radius=r1)
        nodes, indexes, distances = cm.n_closest_points([25, 35, 41], 3)
        self.assertEqual(False, True)

    def test_bbox_collision(self):
        bb1 = [[10, 20], [10, 20], [10, 20]]
        bb2 = [[15, 30], [15, 30], [15, 30]]
        bb3 = [[15, 30], [15, 30], [25, 30]]
        bb4 = [[5, 35], [5, 35], [15, 18]]
        params = [
            [bb1, bb2, True],
            [bb1, bb3, False],
            [bb1, bb4, True],
        ]
        for param in params:
            # out = g3.check_Collision(param[0], param[1])
            out = g3.bbox_collision(param[0], param[1])
            self.assertEqual(out, param[2])

        # g3.bbox_collision(bb1, bb2)

    def test_get_bbox(self):
        points = [[10, 20, 0], [30, 40, 0]]
        bbox_expected = np.asarray([[5, 35], [15, 45], [-5, 5]])

        bbox = g3.get_bbox(points, margin=5)
        err = np.sum((bbox - bbox_expected)**2)
        self.assertAlmostEquals(err, 0)


    def test_dist_between_lines(self):
        # TODO
        a0 = [0, 0, 0]
        a1 = [0, 0, 1]
        b0 = [5, 0, 1]
        b1 = [8, -3, -5]

        pa1, pb1, dist1 = g3.closest_distance_between_lines(a0, a1, b0, b1)
        pa2, pb2, dist2 = g3.closest_distance_between_lines(a0, a1, b0, b1, clampAll=True)

        # self.assertAlmostEquals(dist, 0

    def test_dist_between_lines_problematic_case(self):
        # TODO
        a0 = [15.66021867, 18.12350593, 14.89322993]
        a1 = [15.66021867, 17.10675439, 13.87647838]
        b0 = [16.05191401, 21.03424487, 10.59026252]
        b1 = [16.05191401, 13.6576145, 3.21363215]

        pa1, pb1, dist1 = g3.closest_distance_between_lines(a0, a1, b0, b1)
        pa2, pb2, dist2 = g3.closest_distance_between_lines(a0, a1, b0, b1, clampAll=True)

    def test_dist_between_paralel_lines(self):
        a0 = [0, 0, 0]
        a1 = [0, 0, 1]
        b0 = [1, 0, 0]
        b1 = [1, 0, 1]

        pa, pb, dist = g3.closest_distance_between_lines(a0, a1, b0, b1)

        self.assertAlmostEquals(dist, 0)


    def test_point_in_plane(self):
        pl_pt = [0, 0, 0]
        pl_vect = [1, 0, 0]

        pts = np.asarray([
            [0, 0, 1],
            [0, -1, 1],
            [1, 1, 1],
            [1, 15, 13],
            [-1, -1, 1],
            [-8, -1, 1]
            ]).T

        retval_expected = [
            0,
            0,
            1,
            1,
            -1,
            -1,
        ]
        retval = g3.point_and_plane(pl_pt, pl_vect, pts)
        retval = np.sign(retval)
        err = np.sum((retval - retval_expected)**2)

        self.assertAlmostEqual(err, 0)

    def test_cylinder_collision_model_collision(self):
        cm = g3.CollisionModelSpheres(areasize=[150, 151, 155])
        distance = 25
        pt1 = [20, 20, 20]
        pt2 = [20, 20, 60]
        pt3 = [20, 20 + distance, 20]
        pt4 = [20, 20 + distance, 60]
        ptC1 = [0, 0, 0]
        ptC2 = [0, 60, 0]
        rA = 10
        rB = 10
        rC = 15
        cylA = g3.CylinderObject(pt1, pt2, radius=rA)
        cylB = g3.CylinderObject(pt3, pt4, radius=rB)
        cylC = g3.CylinderObject(ptC1, ptC2, radius=rC)
        collision1 = cylA.collision(cylB)
        collision2 = cylA.collision(cylC)

        self.assertEqual(collision1, False)
        self.assertEqual(collision2, False)

    def test_bbox_corners(self):
        bbox = [[5, 10], [13, 17], [16, 18]]
        points = g3.get_bbox_corners(bbox)

        expected_point = np.array([10,13,16])
        found = False
        for point in points:
            err = np.sum((point - expected_point)**2)
            if err == 0:
                found = True
                break

        self.assertTrue(found)

    def test_cylinder_collision_model_out_of_area(self):
        """
        Cylinder end point should be in safe distance from boundary
        :return:
        """
        cm = g3.CollisionModelSpheres(areasize=[150, 151, 65])
        distance = 25
        pt1 = [20, 20, 20]
        pt2 = [20, 20, 60]
        pt3 = [20, 20 + distance, 20]
        pt4 = [20, 20 + distance, 60]
        r1 = 10
        r2 = 10
        collision1 = cm.add_cylinder_if_no_collision(pt1, pt2, radius=r1)
        collision2 = cm.add_cylinder_if_no_collision(pt3, pt4, radius=r2)

        self.assertEqual(collision2, True)

    def test_cylinder_collision_model_dist_problem(self):
        """
        check where is problem in two lines dist function
        :return:
        """

        cm = g3.CollisionModelSpheres(areasize=[150, 151, 155])

        distance = 25
        pt1 = [20, 20, 20]
        pt2 = [20, 20, 60]
        pt3 = [20, 20 + distance, 20]
        pt4 = [20, 20 + distance, 60]
        ptC1 = [0, 0, 0]
        ptC2 = [0, 60, 0]
        rA = 10
        rB = 10
        rC = 15
        cylA = g3.CylinderObject(
            [ 8.2932388 , -0.65185999,  9.54170155],
            [8.2932388 ,  11.23408448,  21.42764602],
            radius=5.691834898023455)
        cylB = g3.CylinderObject(
            [ 2.76365903,  0.86076129,  4.30361747],
            [2.76365903,   7.00253318,  10.44538936],
            radius=7.2248666037231875
        )
        collision1 = cylA.collision(cylB)

        self.assertEqual(collision1, True)

if __name__ == '__main__':
    unittest.main()
