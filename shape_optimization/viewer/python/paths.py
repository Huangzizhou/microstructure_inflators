#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

# Code adapted from quadfoam project, authored by Jeremie Dumas

import os

MICRO_DIR = "/Users/davitozoni/Desktop/NYU/Research/Repositories/microstructures"
MESHFEM_DIR = "/Users/davitozoni/Desktop/NYU/Research/Repositories/MeshFEM"
PROJECT_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
CONFIG_DIR = os.path.join(PROJECT_DIR, "config")


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    """
    Mimics the behavior of Unix's 'which' command.

    Parameters
    ----------
    program : string
        Absolute or relative name of the program to find.
    """
    fpath, _fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def find(query, folder):
    """
    Finds a file in the given folder hierarchy. In case multiple files are found,
    returns the most recent one.

    Parameters
    ----------
    query: string
        Name of the file to search for.
    folder: string
        Name of the folder to explore.
    """
    candidates = []
    for subdir, _dirs, files in os.walk(folder):
        for file in files:
            fullpath = os.path.join(os.path.abspath(subdir), file)
            basename = os.path.splitext(file)[0]
            if basename == query and is_exe(fullpath):
                candidates.append(fullpath)
    candidates.sort(key=os.path.getmtime)
    if not candidates:
        return None
    return candidates[-1]


def get_executable(query, folders=None):
    """
    Retrieve a program path, either in system paths (via 'which') or in a given
    list of folders (via 'find'). Raise a FileNotFoundError exception in case of
    failure.
    """
    prog = which(query)
    if folders is None:
        folders = []
    if isinstance(folders, str):
        folders = [folders]
    for folder in folders:
        if prog is None:
            prog = find(query, folder)
    if prog is None:
        raise FileNotFoundError(query)
    return prog
