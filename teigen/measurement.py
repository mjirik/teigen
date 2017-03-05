#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

import logging

logger = logging.getLogger(__name__)

import numpy as np
import skimage.measure

def surface_measurement(volume, voxelsize, level=1.0, **kwargs):
    vertices, faces = skimage.measure.marching_cubes(volume, level=level, spacing=voxelsize)
    surface_area = skimage.measure.mesh_surface_area(verts=vertices, faces=faces)
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh = Poly3DCollection(vertices[faces])
    mesh.set_edgecolor('k')
    ax.add_collection3d(mesh)
    sh = volume.shape
    ax.set_xlim(0, sh[0])  # a = 6 (times two for 2nd ellipsoid)
    ax.set_ylim(0, sh[1])  # b = 10
    ax.set_zlim(0, sh[2])  # c = 16
    plt.show()
    return surface_area


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
