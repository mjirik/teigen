#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8


import numpy as np
import skimage.io
import os
import os.path as op


def save_image_stack(data3d, filename="output{:05d}.png"):
    """

    :param data3d:
    :param filename:
    :return:
    """

    for i in range(data3d.shape[0]):
        fn = filename.format(i)
        dn = op.dirname(fn)
        if not op.exists(dn):
            os.makedirs(dn)

        slice = data3d[i]
        # import ipdb; ipdb.set_trace()

        skimage.io.imsave(fn, slice)
