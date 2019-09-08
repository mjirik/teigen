#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © %YEAR%  <>
#
# Distributed under terms of the %LICENSE% license.



"""



"""



import logging




logger = logging.getLogger(__name__)

import argparse


import PyQt5
from PyQt5.QtWidgets import (QGridLayout, QLabel, QPushButton, QApplication,
                         QWidget, QGridLayout, QSpinBox, QLineEdit, QCheckBox,
                         QComboBox, QTextEdit, QDialog, QMainWindow,
                         QDoubleSpinBox)


from PyQt5 import QtGui, QtWidgets
import sys
import os.path as op
import copy
import numpy as np




def str_format_old_to_new(string):

    """

    convert old format style to new style. Works for digits only

    %05d is converted to {:05d}

    :param string:

    :return:

    """

    import re


    return re.sub(r"%(\d*d)", r"{:\1}", string)





class SetDirWidget(QtWidgets.QWidget):

    def __init__(self, input_dir="", caption=None, tooltip=None):

        """



        :param config_in:  dictionary

        :param ncols:

        :param captions:

        """

        super(SetDirWidget, self).__init__()



        self.caption = caption

        self.tooltip = None

        self.input_dir = op.expanduser(input_dir)

        self.init_ui()



    def init_ui(self):

        self.mainLayout = QGridLayout(self)

        self.widgets = {}

        grid = self.mainLayout



        if self.caption is not None:

            grid.addWidget(QLabel(self.caption), 0, 0)



        self.ui_dir_box = QLineEdit()

        self.ui_dir_box.setText(self.input_dir)

        # self.config.add_handler(key, dir_box)

        grid.addWidget(self.ui_dir_box, 0, 1)



        btn_accept = QPushButton("Set dir", self)

        if self.tooltip is not None:

            btn_accept.setToolTip(self.tooltip)

        btn_accept.clicked.connect(self.callback_set_dir)

        grid.addWidget(btn_accept, 0, 2)



    def get_dir(self):

        return op.expanduser(str(self.ui_dir_box.text()))



    def callback_set_dir(self):

        # init_filename = "file%05d.jpg"

        directory, file = op.split(self.get_dir())



        filename = QtWidgets.QFileDialog.getExistingDirectory(

            self,

            "Save directory",

            directory

        )

        filename = op.join(str(filename), file)



        self.ui_dir_box.setText(filename)





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

