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
    if type(dct) == collections.OrderedDict:
        p = collections.OrderedDict()
    else:
        p = {}
    for key, value in dct.items():
        if key in keys:
            p[key] = value
    # p = {key: value for key, value in dct.items() if key in keys}
    return p

def kick_from_dict(dct, keys):
    if type(dct) == collections.OrderedDict:
        p = collections.OrderedDict()
    else:
        p = {}
    for key, value in dct.items():
        if key not in keys:
            p[key] = value

    # p = {key: value for key, value in dct.items() if key not in keys}
    return p

def split_dict(dct, keys):
    """
    Split dict into two subdicts based on keys
    :param dct:
    :param keys:
    :return: dict_in, dict_out
    """
    if type(dct) == collections.OrderedDict:
        dict_in = collections.OrderedDict()
        dict_out = collections.OrderedDict()
    else:
        dict_in = {}
        dict_out = {}

    for key, value in dct.items:
        if key in keys:
            dict_in[key] = value
        else:
            dict_out[key] = value
    return dict_in, dict_out

def recursive_update(d, u):
    """
    Dict recursive update.

    Based on Alex Martelli code on stackoverflow
    http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth?answertab=votes#tab-top

    :param d:
    :param u:
    :return:
    """
    for k, v in u.iteritems():
        if isinstance(v, collections.Mapping):
            r = recursive_update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d