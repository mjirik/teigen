#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

import logging

logger = logging.getLogger(__name__)

import numpy as np
import skimage.measure


def surface_measurement(volume, voxelsize, level=None, return_vertices_and_faces=False, **kwargs):
    if level==None:
        level = (np.max(volume) + np.min(volume)) * 0.5
    vertices, faces = skimage.measure.marching_cubes(volume, level=level, spacing=voxelsize)
    surface_area = skimage.measure.mesh_surface_area(verts=vertices, faces=faces)

    if return_vertices_and_faces:
        return surface_area, vertices, faces
    else:
        return surface_area


def volume_measurement(volume, voxelsize, level=1.0, return_vertices_and_faces=False, **kwargs):
    volume = np.sum(volume > level) * np.prod(voxelsize)
    return volume


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
