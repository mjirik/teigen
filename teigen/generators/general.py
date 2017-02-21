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

class GeneralGenerator:

    def __init__(self):

        self.data3d = None

    def generate_volume(self, *args, **kwargs):
        from ..tree import TreeBuilder

        self.tvgvol = TreeBuilder('vol')
        self.tvgvol.voxelsize_mm = self.voxelsize_mm # [1, 1, 1]
        self.tvgvol.shape = self.areasize_px # [100, 100, 100]
        self.tvgvol.tree_data = self.tree_data
        self.tvgvol.finish_progress_callback = self.progress_callback
        if self.intensity_profile is not None:
            self.tvgvol.intensity_profile = self.intensity_profile
        self.data3d = self.tvgvol.buildTree(*args, **kwargs)
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

