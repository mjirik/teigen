#! /usr/bin/python
# -*- coding: utf-8 -*-

"""

"""
import logging

logger = logging.getLogger(__name__)

# import funkcí z jiného adresáře
import sys
import os.path
import scipy


def filter_data(data3d, voxelsize_mm, gaussian_filter=True, gaussian_filter_sigma_mm=1.0):
    sigma_px = gaussian_filter_sigma_mm / voxelsize_mm
    data3d = scipy.ndimage.filters.gaussian_filter(
        data3d,
        sigma=sigma_px)

    return data3d


def main():
    logger = logging.getLogger()

    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
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
        '-i', '--inputfile',
        default=None,
        required=True,
        help='input file'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='Debug mode')
    args = parser.parse_args()

    if args.debug:
        ch.setLevel(logging.DEBUG)

    app = QApplication(sys.argv)
    # cfg = {"bool": True, "int":5, 'str': 'strdrr'}
    # captions = {"int": "toto je int"}
    cw = SetDirWidget("~/lisa_data")
    cw.show()
    app.exec_()


if __name__ == "__main__":
    main()
