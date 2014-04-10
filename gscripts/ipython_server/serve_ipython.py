#!/usr/bin/python
import sys
from subprocess import Popen

port = sys.argv[1]
mongoport=sys.argv[2]
print port, mongoport

try:
    c = Popen(["ssh", "-L", ("%s:localhost:%s" %(mongoport, mongoport)), "oolite", "-N"])
except:
    print "not linking mongodb"


b = Popen(["ipython2", "notebook", "--no-browser", "--port", "%s" %(port), "--pylab", "inline"])

a = Popen(["ssh", "-R", ("%s:localhost:%s" %(port, port)), "oolite", "-N"])
b.wait()
a.wait()

if c:
    c.wait()
