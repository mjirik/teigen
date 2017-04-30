#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

import logging

logger = logging.getLogger(__name__)

import argparse

import vtk

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


def main():
    logger = logging.getLogger()

    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    logger.addHandler(ch)

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
        default=None,
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


if __name__ == "__main__":
    main()
