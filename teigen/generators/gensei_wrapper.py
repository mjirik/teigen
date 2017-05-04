#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

import logging

logger = logging.getLogger(__name__)

import os, glob, copy, time
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter


class GenseiGenerator:
    def __init__(self,
                 # voxelsize_mm=[1.0, 1.0, 1.0],
                 # area_shape=[100, 100, 100],
                 resolution=[600, 600],
                 dims=[10, 10, 10],
                 n_slice=21,
                 n_objects=10
                 ):
        self.objects = None
        # self.voxelsize_mm = np.asarray(voxelsize_mm)
        # self.area_shape = np.asarray(area_shape)
        # self.area_volume = np.prod(self.area_shape * self.voxelsize_mm)
        self.resolution = np.asarray(resolution)
        self.dims = dims = np.asarray(dims)
        self.n_slice = n_slice
        self.n_objects = n_objects
        self.area_volume = np.prod(self.dims)

        pass

    def run(self):
        logger.info("gensei generator running")
        from gensei.base import Config, Object
        from gensei import Objects, Box
        from gensei.utils import get_suffix

        conf = {'objects': default_objects,
                'box': default_box,
                'options': default_options}

        config = Config.from_conf(conf, required=['objects', 'box'],
                                  optional=['options'])

        if isinstance(config.box['dims'], str):
            config.box['dims'] = eval(config.box['dims'])
        if isinstance(config.box['resolution'], str):
            aux = tuple([int(r) for r in config.box['resolution'].split('x')])
            config.box['resolution'] = aux

        # config.box["resolution"] = 1/self.voxelsize_mm[1:]
        # config.box['dims'] = self.area_shape
        config.box["resolution"] = tuple(self.resolution.astype(int).tolist())
        config.box['dims'] = self.dims.astype(int).tolist()
        config.box['n_slice'] = int(self.n_slice)
        box = Box(**config.box)
        options = Object(name='options', **config.options)

        output(box)
        output(options)

        object_classes = Objects.from_conf(config.objects, box)
        objects = object_classes.place_objects(box, options)

        self.objects = objects
        self.box = box
        self.options = options
        pass

    def get_stats(self):
        # self.objects.format_intersection_statistics(is_output=True)
        # output('saving report to %s' % reportname)
        # fd = open(reportname, 'w')
        #
        # fd.write('started: %s\n' % time.ctime(time_start))
        # fd.write('elapsed: %.1f [s]\n' % (time_end - time_start))
        #
        # box.report(fd)
        # options.report(fd)
        #
        # fd.write('\n'.join(objects.format_statistics())+'\n')
        #
        # object_classes.report(fd)
        # objects.report(fd)
        #
        # fd.close()
        #
        # output('all done.')
        data_dict = {
            'volume': [],
            'surface': [],
            'length': [],
            'radius': [],
        }
        for key in self.objects:
            obj = self.objects[key]
            data_dict['surface'].append(obj.surface)
            data_dict['volume'].append(obj.volume)
            data_dict['length'].append(0)
            data_dict['radius'].append(0)
            # obj.

        import pandas as pd
        df = pd.DataFrame(data_dict)
        return df

    def saveVolumeToFile(self, dirname):
        if self.objects is None:
            self.run()

        generate_slices(self.objects, self.box, self.options, dirname)


def generate_slices(objects, box, options, output_filename_trunk):
    """
    Save images of the specimen slices along the specified axes of the
    block. Each image displays a planar cut plane of the block intersecting the
    ellipsoids.
    """
    from gensei.base import output, assert_
    from gensei import Objects, Box
    from gensei.utils import get_suffix
    resolution = box.resolution

    imshape = resolution[::-1] + (3,)
    aspect = float(resolution[1]) / resolution[0]
    figsize = plt.figaspect(aspect)
    dpi = resolution[0] / figsize[0]

    objects.init_section_based_data()
    objects.points = []
    for pb, points, delta, n_slice, axis, am in box.get_points():
        suffix = get_suffix(n_slice)

        # dpi=dpi in plt.figure() messes with figsize... ???
        fig = plt.figure(1, figsize=figsize, dpi=dpi)
        fig.set_figwidth(figsize[0])
        fig.set_figheight(figsize[1])
        ax = fig.add_axes([0, 0, 1, 1])

        objects.init_section_based_data(axis)

        x1b, x2b = pb[am[0]], pb[am[1]]
        for islice, x3b in enumerate(pb[am[2]]):
            x3b_name = ('%05.2f' % x3b).replace('.', '_')
            filename = '.'.join((output_filename_trunk, axis,
                                 suffix % islice, x3b_name,
                                 options.output_format))
            output(islice, x3b, filename, '...')
            output('computing')
            points[:, am[2]] = x3b

            mask = np.zeros(points.shape[0], dtype=np.int8)
            cmask = np.zeros((points.shape[0], 3), dtype=np.float64)
            for obj in objects.itervalues():
                color = np.array(colorConverter.to_rgb(obj.conf.color))

                bbox = obj.get_aligned_bounding_box()[am]

                ix1 = np.where((x1b > bbox[0, 0]) & (x1b < bbox[0, 1]))[0]
                ix2 = np.where((x2b > bbox[1, 0]) & (x2b < bbox[1, 1]))[0]
                a, b = np.meshgrid(ix1, resolution[0] * ix2)
                ii = (a + b).ravel()

                _mask = obj.contains(points[ii])
                mask[ii] += _mask
                cmask[ii[_mask]] = color
                objects.update_section_based_data(_mask, a.shape, axis, delta,
                                                  islice, x3b, obj)
                obj.store_intersection(_mask, axis, x3b)

            objects.points.append((axis, islice, x3b))

            try:
                assert_(np.alltrue(mask <= 1))
            except:
                import pdb;
                pdb.set_trace()

            output('drawing')
            ax.cla()
            ax.set_axis_off()
            ax.imshow(cmask.reshape(imshape), origin='upper')

            output('saving')
            plt.savefig(filename, format=options.output_format, dpi=dpi)
            output('...done')


##        plt.show()



objects = {
    'class 1': {
        'kind': 'ellipsoid',
        'color': 'r',
        'fraction': 0.1,
        'length_to_width': 2.0,
        'reduce_to_fit': {'length_to_width': 0.9},
        'centre': 'random',
        # Ellipsoids should be approximately aligned with the x axis.
        'rot_axis': [('normal', 0.0, 0.1), 1.0, ('random', -0.2, 0.2)],
        'rot_angle': np.pi / 2,
    },
    'class 2': {
        'kind': 'cylinder',
        'color': (0.1, 0.2, 0.7),
        'fraction': 0.01,
        'length_to_width': 40.0,
        'reduce_to_fit': {'fraction': 0.9},
        'centre': 'random',
        # Cylinders should have uniform random distribution of directions.
        'direction': 'random direction',
    },
}
# --------End of settings of the properties of objects-----------

# --------Settings of the properties of the box-----------
box = {
    # dimensions of the box
    'dims': (10, 10, 10),
    # arbitrary units
    'units': 'mm',
    # size of the generated image in pixels
    'resolution': (300, 300),
    # number of objects to be generated within the box, either a number, or a
    # per class dictionary, for example:
    'n_object': {'class 1': 10, 'class 2': 50},
    # number of slices generated - a dictionary of numbers for the
    # directions perpendicular to the Z, X, and Y axis, i.e. in the XY, YZ, and
    # XY planes, or a single integer.
    'n_slice': {'z': 21, 'x': 21, 'y': 21},
}
# --------End of settings of the properties of objects-----------

# --------General options-----------
options = {
    # output file format supported by the matplotlib backend used
    'output_format': 'png',
    # timeout in seconds to place more ellipsoids in to the box
    'timeout': 5.0,
}

defaults = {
    'fraction': 0.1,
    'fraction_reduction': 0.9,
    'length_to_width': 8.0,

    'n_slice': 21,
    'dims': '(10, 10, 10)',
    'units': 'mm',
    'resolution': '600x600',
    'n_object': 10,
    'output_format': 'png',
    'timeout': 5.0,
}

default_objects = {
    'class 1': {
        'kind': 'ellipsoid',
        'color': 'r',
        'fraction': defaults['fraction'],
        'length_to_width': defaults['length_to_width'],
        'reduce_to_fit': {'fraction': defaults['fraction_reduction']},
        'centre': 'random',
        'direction': 'random direction',
    },
}

default_box = {
    'dims': defaults['dims'],
    'units': defaults['units'],
    'resolution': defaults['resolution'],
    'n_object': defaults['n_object'],
    'n_slice': defaults['n_slice'],
}

default_options = {
    'output_format': defaults['output_format'],
    'timeout': defaults['timeout'],
}
