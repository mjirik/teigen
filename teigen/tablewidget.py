# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 10:50:56 2013
@author: cmarshall
"""

import logging
logger = logging.getLogger(__name__)
import sip

sip.setapi('QString', 1)
sip.setapi('QVariant', 1)

import pandas as pd
from PyQt5 import QtCore, QtGui, QtWidgets


class TableModel(QtCore.QAbstractTableModel):
    def __init__(self, parent=None, *args):
        super(TableModel, self).__init__()
        self.datatable = None

    def update(self, dataIn):
        logger.info('Updating Model')
        self.datatable = dataIn
        self.headerdata = list(dataIn.keys())
        self.indexdata = list(dataIn.index)

        logger.debug('Datatable : {0}'.format(self.datatable))
        logger.debug(self.headerdata)
        logger.debug(self.indexdata)



    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self.datatable.index)

    def columnCount(self, parent=QtCore.QModelIndex()):

        return len(self.datatable.columns.values)



    def data(self, index, role=QtCore.Qt.DisplayRole):
        # logger.debug('Data Call')
        # print index.column(), index.row()

        if role == QtCore.Qt.DisplayRole:
            i = index.row()
            j = index.column()
            # return QtCore.QVariant(str(self.datatable.iget_value(i, j)))
            # return '{0}'.format(self.datatable.iget_value(i, j))
            return '{0}'.format(self.datatable.iat[i, j])
        else:
            return QtCore.QVariant()

    def flags(self, index):
        return QtCore.Qt.ItemIsEnabled

    def headerData(self, col, orientation, role):
        # return str(self.headerdata[col])
        from PyQt5.QtCore import QVariant, Qt

        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headerdata[col])
        return QVariant()

class TableView(QtWidgets.QTableView):
    """
    A simple table to demonstrate the QComboBox delegate.
    """
    def __init__(self, *args, **kwargs):
        QtWidgets.QTableView.__init__(self, *args, **kwargs)


class TableWidget(QtWidgets.QWidget):
    """
    A simple test widget to contain and own the model and table.
    """



    def __init__(self, parent=None, dataframe=None):

        QtWidgets.QWidget.__init__(self, parent)
        l = QtWidgets.QVBoxLayout(self)
        if dataframe is None:
            cdf = self.get_data_frame()
        else:
            cdf = dataframe
        self._tm = TableModel(self)
        self._tm.update(cdf)
        self._tv = TableView(self)
        self._tv.setModel(self._tm)
        l.addWidget(self._tv)

    def get_data_frame(self):
        df = pd.DataFrame({'Name': ['a', 'b', 'c', 'd'],
                           'First': [2.3, 5.4, 3.1, 7.7], 'Last': [23.4, 11.2, 65.3, 88.8], 'Class': [1, 1, 2, 1],
                           'Valid': [True, True, True, False]})
        return df


if __name__ == "__main__":
    from sys import argv, exit
    a = QtWidgets.QApplication(argv)
    w = TableWidget()
    w.show()
    w.raise_()
    exit(a.exec_())
