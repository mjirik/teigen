# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)
import unittest
import teigen.geometry3d as g3

class GeometryTestCase(unittest.TestCase):
    def test_collision_model_no_collision(self):
        cm = g3.CollisionBoundaryModel(areasize=[150, 151, 155])
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
        cm = g3.CollisionBoundaryModel(areasize=[150, 151, 155])
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
        cm = g3.CollisionBoundaryModel(areasize=[150, 151, 65])
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

if __name__ == '__main__':
    unittest.main()
