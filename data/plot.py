"""
plot norms v. time
"""

# imports
import matplotlib.pyplot as plt
import sys
import getopt

#opts, dumps = getopt.getopt(sys.argv[1:], "-f:-i:")
#for opt, arg in opts:
#    if opt == "-f":
#        filename = str(arg);
#    if opt == "-i":
#        value_idx = int(arg);
 
filename = sys.argv[1]
value_idx = int(sys.argv[2])
time_list = []
value_list = []

with open(filename,"r") as f:
 for line in f.readlines():
     obj_list = line.split(" ")
     time_list.append(float(obj_list[0]))
     value_list.append(float(obj_list[value_idx]))

plt.xscale('log')
plt.yscale('log')
plt.plot(time_list, value_list, 'r')

plt.show()
