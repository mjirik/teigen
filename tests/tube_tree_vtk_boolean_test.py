#! /usr/bin/python
# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)
import unittest
import pytest

import os
import os.path

import vtk
import teigen.tb_vtk
import numpy as np


class VtkBooleanTestCase(unittest.TestCase):
    @pytest.mark.interactive
    def test_vtk_example_with_my_cylinder(self):
        # sphereSource1 = vtk.vtkSphereSource()
        # sphereSource1.SetCenter(0.25, 0, 0)
        # sphereSource1.Update()
        # input1 = sphereSource1.GetOutput()
        input1 = teigen.tb_vtk.get_cylinder([0.25, 0, -.5], 0.9, 0.7, [0.0, .0, .0])
        sphere1tri = vtk.vtkTriangleFilter()
        sphere1tri.SetInputData(input1)

        # sphereSource2 = vtk.vtkSphereSource()
        # sphereSource2 = vtk.vtkCylinderSource()
        # sphereSource2.Update()
        # input2 = sphereSource2.GetOutput()
        direction = np.asarray([1., 1., 1.])
        direction /= np.linalg.norm(direction)
        input2 = teigen.tb_vtk.get_cylinder([0, .1, 0], 0.7, 0.5, direction)

        sphere2Tri = vtk.vtkTriangleFilter()
        sphere2Tri.SetInputData(input2)

        input1mapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            input1mapper.SetInputConnection(input1.GetProducerPort())
        else:
            input1mapper.SetInputData(input1)

        input1mapper.ScalarVisibilityOff()
        input1Actor = vtk.vtkActor()
        input1Actor.SetMapper(input1mapper)
        input1Actor.GetProperty().SetColor(1, 0, 0)
        input1Actor.SetPosition(input1.GetBounds()[1] - input1.GetBounds()[0], 0, 0)
        input2Mapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            input2Mapper.SetInputConnection(input2.GetProducerPort())
        else:
            input2Mapper.SetInputData(input2)

        input2Mapper.ScalarVisibilityOff()
        input2Actor = vtk.vtkActor()
        input2Actor.SetMapper(input2Mapper)
        input2Actor.GetProperty().SetColor(0, 1, 0)
        input2Actor.SetPosition(-(input2.GetBounds()[1] - input2.GetBounds()[0]), 0, 0)

        booleanOperation = vtk.vtkBooleanOperationPolyDataFilter()
        # booleanOperation.SetOperationToUnion()
        # booleanOperation.SetOperationToIntersection()
        booleanOperation.SetOperationToDifference()

        if vtk.VTK_MAJOR_VERSION <= 5:
            booleanOperation.SetInputConnection(0, sphere1tri.GetOutputPort())
            booleanOperation.SetInputConnection(1, sphere2Tri.GetOutputPort())
        else:
            sphere1tri.Update()
            sphere2Tri.Update()
            booleanOperation.SetInputData(0, sphere1tri.GetOutput())
            booleanOperation.SetInputData(1, sphere2Tri.GetOutput())
        booleanOperation.Update()

        boolean_peration_mapper = vtk.vtkPolyDataMapper()
        boolean_peration_mapper.SetInputConnection(booleanOperation.GetOutputPort())
        boolean_peration_mapper.ScalarVisibilityOff()

        booleanOperationActor = vtk.vtkActor()
        booleanOperationActor.SetMapper(boolean_peration_mapper)

        renderer = vtk.vtkRenderer()
        renderer.AddViewProp(input1Actor)
        renderer.AddViewProp(input2Actor)
        renderer.AddViewProp(booleanOperationActor)
        renderer.SetBackground(.1, .2, .3)
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.AddRenderer(renderer)

        ren_win_interactor = vtk.vtkRenderWindowInteractor()
        ren_win_interactor.SetRenderWindow(renderWindow)

        renderWindow.Render()
        ren_win_interactor.Start()
        # self.assertEqual(True, False)

    @pytest.mark.interactive
    def test_vtk_example(self):
        sphereSource1 = vtk.vtkSphereSource()
        sphereSource1.SetCenter(0.25, 0, 0)
        sphereSource1.Update()
        input1 = sphereSource1.GetOutput()
        sphere1Tri = vtk.vtkTriangleFilter()
        sphere1Tri.SetInputData(input1)

        # sphereSource2 = vtk.vtkSphereSource()
        sphereSource2 = vtk.vtkCylinderSource()
        sphereSource2.Update()
        input2 = sphereSource2.GetOutput()
        # input2 = teigen.tb_vtk.get_cylinder([0, 0, 0], 0.5, 0.5, [1,1,1])
        sphere2Tri = vtk.vtkTriangleFilter()
        sphere2Tri.SetInputData(input2)

        input1Mapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            input1Mapper.SetInputConnection(input1.GetProducerPort())
        else:
            input1Mapper.SetInputData(input1)

        input1Mapper.ScalarVisibilityOff()
        input1Actor = vtk.vtkActor()
        input1Actor.SetMapper(input1Mapper)
        input1Actor.GetProperty().SetColor(1, 0, 0)
        input1Actor.SetPosition(input1.GetBounds()[1] - input1.GetBounds()[0], 0, 0)
        input2Mapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            input2Mapper.SetInputConnection(input2.GetProducerPort())
        else:
            input2Mapper.SetInputData(input2)

        input2Mapper.ScalarVisibilityOff()
        input2Actor = vtk.vtkActor()
        input2Actor.SetMapper(input2Mapper)
        input2Actor.GetProperty().SetColor(0, 1, 0)
        input2Actor.SetPosition(-(input2.GetBounds()[1] - input2.GetBounds()[0]), 0, 0)

        booleanOperation = vtk.vtkBooleanOperationPolyDataFilter()
        # booleanOperation.SetOperationToUnion()
        # booleanOperation.SetOperationToIntersection()
        booleanOperation.SetOperationToDifference()

        if vtk.VTK_MAJOR_VERSION <= 5:
            booleanOperation.SetInputConnection(0, sphere1Tri.GetOutputPort())
            booleanOperation.SetInputConnection(1, sphere2Tri.GetOutputPort())
        else:
            sphere1Tri.Update()
            sphere2Tri.Update()
            booleanOperation.SetInputData(0, sphere1Tri.GetOutput())
            booleanOperation.SetInputData(1, sphere2Tri.GetOutput())
        booleanOperation.Update()

        booleanOperationMapper = vtk.vtkPolyDataMapper()
        booleanOperationMapper.SetInputConnection(booleanOperation.GetOutputPort())
        booleanOperationMapper.ScalarVisibilityOff()

        booleanOperationActor = vtk.vtkActor()
        booleanOperationActor.SetMapper(booleanOperationMapper)

        renderer = vtk.vtkRenderer()
        renderer.AddViewProp(input1Actor)
        renderer.AddViewProp(input2Actor)
        renderer.AddViewProp(booleanOperationActor)
        renderer.SetBackground(.1, .2, .3)
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.AddRenderer(renderer)

        renWinInteractor = vtk.vtkRenderWindowInteractor()
        renWinInteractor.SetRenderWindow(renderWindow)

        renderWindow.Render()
        renWinInteractor.Start()
        # self.assertEqual(True, False)

    # @pytest.mark.interactive
    def test_vtk_surface_and_volume(self):
        import teigen.geometry3d as g3
        height = 1.0
        radius = 1.0
        input1 = teigen.tb_vtk.get_cylinder([0.25, 0, -.5],
                                            height=height,
                                            radius=radius,
                                            direction=[0.0, .0, .0],
                                            resolution=50
                                            )
        object1Tri = vtk.vtkTriangleFilter()
        object1Tri.SetInputData(input1)
        object1Tri.Update()
        mass = vtk.vtkMassProperties()
        mass.SetInputData(object1Tri.GetOutput())
        surf = mass.GetSurfaceArea()
        vol = mass.GetVolume()

        surf_analytic = g3.cylinder_surface(radius, height)
        vol_analytic = g3.cylinder_volume(radius, height)
        err_surf = np.abs(surf_analytic - surf) / surf_analytic
        err_vol = np.abs(vol_analytic - vol) / vol_analytic
        # logger.debug(surf, surf_analytic, err_surf)
        # logger.debug(vol, vol_analytic, err_vol)

        max_error = 0.01
        self.assertLess(err_surf, max_error)
        self.assertLess(err_vol, max_error)

    def test_vtk_check_orientation(self):
        import teigen.geometry3d as g3
        pt1 = np.asarray([4.22, 8.56, 9.39])
        pt2 = np.asarray([3.24, 9.46, 2.83])
        height = np.linalg.norm(pt1 - pt2)

        direction=(pt1 - pt2) / height

        radius = 2.53
        input1 = teigen.tb_vtk.get_cylinder(pt1,
                                            height=height,
                                            radius=radius,
                                            direction=direction,
                                            resolution=25
                                            )
        bounds1 = input1.GetBounds()
        object1Tri = vtk.vtkTriangleFilter()
        object1Tri.SetInputData(input1)
        object1Tri.Update()
        mass = vtk.vtkMassProperties()
        mass.SetInputData(object1Tri.GetOutput())
        bounds0 = object1Tri.GetOutput().GetBounds()

        surf = mass.GetSurfaceArea()
        vol = mass.GetVolume()

        surf_analytic = g3.cylinder_surface(radius, height)
        vol_analytic = g3.cylinder_volume(radius, height)
        err_surf = np.abs(surf_analytic - surf) / surf_analytic
        err_vol = np.abs(vol_analytic - vol) / vol_analytic
        # logger.debug(surf, surf_analytic, err_surf)
        # logger.debug(vol, vol_analytic, err_vol)

        input2 = teigen.tb_vtk.get_tube(point=pt1,
                                        length=height,
                                        radius=radius,
                                        direction=direction,
                                        sphere_resolution=25,
                                        cylinder_resolution=25
                                        )
        bounds2 = input2.GetBounds()
        max_error = 0.05
        self.assertLess(err_surf, max_error)
        self.assertLess(err_vol, max_error)

if __name__ == '__main__':
    unittest.main()
