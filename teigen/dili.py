#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module is provides funcions for dict lists and functions processing
"""
import logging
logger = logging.getLogger(__name__)
import collections
import inspect


def get_default_args(obj):
    if ("__init__" in dir(obj)):
        if inspect.isfunction(obj.__init__) or inspect.ismethod(obj.__init__):
            argspec = inspect.getargspec(obj.__init__)
        else:
            argspec = inspect.getargspec(obj)
    else:
        argspec = inspect.getargspec(obj)

    args = argspec.args[1:]
    defaults = argspec.defaults
    dc = collections.OrderedDict(zip(args, defaults))
    return dc

def subdict(dct, keys):
    p = {key: value for key, value in dct.items() if key in keys}
    return p
