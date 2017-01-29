# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)
import unittest
import numpy as np
import teigen.geometry3d as g3

class GeometryTestCase(unittest.TestCase):
    def test_collision_model_no_collision(self):
        cm = g3.CollisionSpheresModel(areasize=[150, 151, 155])
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
        cm = g3.CollisionSpheresModel(areasize=[150, 151, 155])
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
        cm = g3.CollisionSpheresModel(areasize=[150, 151, 65])
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
        cm = g3.CollisionSpheresModel(areasize=[150, 151, 155])
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

        # self.assertAlmostEquals(dist, 0)

    def test_dist_between_paralel_lines(self):
        a0 = [0, 0, 0]
        a1 = [0, 0, 1]
        b0 = [1, 0, 0]
        b1 = [1, 0, 1]

        pa, pb, dist = g3.closest_distance_between_lines(a0, a1, b0, b1)

        self.assertAlmostEquals(dist, 0)
if __name__ == '__main__':
    unittest.main()
