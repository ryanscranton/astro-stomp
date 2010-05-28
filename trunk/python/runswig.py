import subprocess
import sys
from sys import stdout

# figure out how swig should be called
swig_addflags='-DWITH_PYTHON'

# also try to add numpy support
try:
    import numpy
    swig_addflags += ' -DWITH_NUMPY'
except:
    sys.stdout.write("Numpy not found, not building with numpy support\n")


# run swig
swig_call = 'swig %s -python -c++ -o stomp_wrap.cxx stomp.i' % swig_addflags
stdout.write("Command: '%s'\n" % swig_call)
pobj = subprocess.Popen(swig_call,shell=True)
stdout_ret, stderr_ret = pobj.communicate()

