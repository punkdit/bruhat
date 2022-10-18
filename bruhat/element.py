#!/usr/bin/env python

from bruhat.argv import argv

if argv.fast:
    print("import bruhat._element")
    from bruhat._element import *

else:
    from bruhat.slow_element import *


