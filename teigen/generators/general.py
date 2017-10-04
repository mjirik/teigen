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

import numpy as np

class GeneralGenerator(object):

    def __init__(self):

        self.data3d = None
        self.stop = False

    def generate_volume(self, *args, **kwargs):
        from ..tb_volume import TBVolume

        self.tvgvol = TBVolume(**kwargs)
        self.tvgvol.voxelsize_mm = self.voxelsize_mm # [1, 1, 1]
        self.tvgvol.shape = self.areasize_px # [100, 100, 100]
        self.tvgvol.tube_skeleton = self.tree_data
        self.tvgvol.finish_progress_callback = self.progress_callback
        if self.intensity_profile is not None:
            self.tvgvol.intensity_profile = self.intensity_profile
        self.data3d = self.tvgvol.buildTree()
        # self.data3d = self.tvgvol.buildTree(*args, **kwargs)
        return self.data3d


    def saveVolumeToFile(self, filename="output{:06d}.jpg"):
        if self.data3d is None:
            self.generate_volume()

        # self.tvgvol.saveToFile(filename)
        import io3d
        import io3d.misc
        import numpy as np
        data = {
            'data3d': self.data3d.astype(np.uint8), #* self.output_intensity,
            'voxelsize_mm': self.voxelsize_mm,
            # 'segmentation': np.zeros_like(self.data3d, dtype=np.int8)
        }
        io3d.write(data, filename)

    def stop(self):
        self.stop = True


def random_normal(mean, scale, min=.0, **kwargs):
    if scale > 0:
        out = np.random.normal(mean, scale, **kwargs)
    else:
        out = mean
    if out < min:
        out = min

    return out
