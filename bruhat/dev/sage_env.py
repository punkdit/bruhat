#!/usr/bin/env python

from __future__ import print_function

import os
os.environ["SAGE_LOCAL"] = "/usr"
os.environ["SAGE_ROOT"] = '/usr/share/sagemath'
os.environ["SAGE_SRC"] = '/usr/share/sagemath/src'

for (k,v) in {
 'CONWAY_POLYNOMIALS_DATA_DIR': '/usr/share/sagemath/conway_polynomials',
 'DOT_SAGE': '/extra/simon/.sage',
 'ELLCURVE_DATA_DIR': '/usr/share/sagemath/ellcurves',
 'GAP_ROOT_DIR': '/usr/share/gap',
 'GRAPHS_DATA_DIR': '/usr/share/sagemath/graphs',
 'HOSTNAME': 'monic',
 'LOCAL_IDENTIFIER': 'monic.10722',
 'POLYTOPE_DATA_DIR': '/usr/share/sagemath/reflexive_polytopes',
 'PYTHON_EGG_CACHE': '/extra/simon/.sage/.python-eggs',
 'REALM': 'sage.math.washington.edu',
 'SAGE_BANNER': '',
 'SAGE_DATE': '2017-12-07',
 'SAGE_DISTFILES': '/usr/share/sagemath/upstream',
 'SAGE_DOC': '/usr/share/doc/sagemath',
 'SAGE_DOC_SRC': '/usr/share/sagemath/src/doc',
 'SAGE_DOT_GIT': '/usr/share/sagemath/.git',
 'SAGE_ETC': '/usr/etc',
 'SAGE_EXTCODE': '/usr/share/sagemath/ext',
 'SAGE_IMPORTALL': 'yes',
 'SAGE_INC': '/usr/include',
 'SAGE_LIB': '/usr/lib/python2.7/dist-packages',
 'SAGE_LOCAL': '/usr',
 'SAGE_LOGS': '/usr/share/sagemath/logs/pkgs',
 'SAGE_PKGS': '/usr/share/sagemath/build/pkgs',
 'SAGE_REPO_ANONYMOUS': 'git://trac.sagemath.org/sage.git',
 'SAGE_REPO_AUTHENTICATED': 'ssh://git@trac.sagemath.org:2222/sage.git',
 'SAGE_ROOT': '/usr/share/sagemath',
 'SAGE_SCRIPTS_DIR': '/usr/share/sagemath/bin',
 'SAGE_SHARE': '/usr/share/sagemath',
 'SAGE_SPKG_INST': '/usr/share/sagemath/installed',
 'SAGE_SRC': '/usr/share/sagemath/src',
 'SAGE_STARTUP_FILE': '/extra/simon/.sage/init.sage',
 'SAGE_URL': 'http://sage.math.washington.edu/sage/',
 'SAGE_VERSION': '8.1',
 'SINGULAR_SO': '/usr/lib/x86_64-linux-gnu/libsingular-Singular-4.1.0.so',
 #'SITE_PACKAGES': ['/usr/lib/python2.7/dist-packages'],
 'THEBE_DIR': '/usr/share/thebe',
 'TRAC_SERVER_URI': 'https://trac.sagemath.org',
 'UNAME': 'Linux'}.items():
    os.environ[k] = v


from sage.all_cmdline import *   # import sage library

