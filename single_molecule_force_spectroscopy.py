#!/usr/bin/env python
# script written for Python 2.7.3 (default, Nov 26 2012, 23:15:21)

###############################################################################
#             Single molecule pulling simulation written by                   #
#                     Alessandro Pirrotta                                     #
###############################################################################

from sys import argv
import datetime
import os
import numpy as np
import math


def mag(vec):
    return math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)


# MD simulation run with tinker 6
# 'dynamic' 'inputmolecule.xyz' INT fs[1.0] ps[0.1] 2 Kelvin
# 'dynamic' 'inputmolecule.xyz' INT fs[1.0] ps[0.1] 1 kcal/mol
#
# INT: 'Number of Dynamics Steps'
# 'Time Step Length in Femtoseconds'
# 'Time between Dumps in Picoseconds'
# '(1) Constant Total Energy Value (E) or
# (2) Constant Temperature via Thermostat (T)'
# 'energy or temperature in Kelvin'

cwd = os.getcwd()

dl = float(1.0e-1)                # stretch per step (usually 0.01 AA)
Dl = float(2.0)                   # total distance of pulling
n = int(Dl/dl)                    # number of cycle
initial_dist = float((2 - 1)*Dl)  # corresponds to dist. in AA between endatoms
int_time = float(1.0)             # in fs integration time (not larger than 1fs)
# thermalization_time_ns = 5.0
thermalization_time_ns = 0.002
steps_thermalization = int(1000000*thermalization_time_ns/int_time)
vel = float(1.0e-2)               # usually 0.001 AA/ps
dump = float(5.0)                 # dump time between dumped snapshots in fs
ndumpstr = int(dl/(dump*vel))
dumps_time_therm = float(1000*thermalization_time_ns)
mode = 2                          # 1 = NVE; 2 = NVT
temp_ener = 300                   # if mode = 1 then K; if mode = 2 then kcal

# no. of md steps between displacement steps
steps = int(((dump*1000.0*ndumpstr)/int_time))

print """
dl, steps, dump, ndumpstr, int_time, dumps_time_therm, vel, steps_thermalization
"""
print dl, steps, dump, ndumpstr, int_time,
print dumps_time_therm, vel, steps_thermalization

# cantilever force constant (in kcal/A^2) (1kcal/A^2mol = 69.5 pN/A=0.695N/m)
# This is   1.1   N/m = 110 pN/AA
kcantilever = float(1.58273)

thermalization = bool(True)  # always run thermalization before production run
outputdat = cwd+'/output.dat'
out = open(outputdat, 'w')
out.write("""
###############################################################################
####             Single molecule pulling simulation written by             ####
####                     Alessandro Pirrotta                               ####
###############################################################################
""")

# check the number of structrures that will be generated
out.write('Total number of snapshot recorded from MD is %s after thermalization'
          % ((n*2))+'\n')
out.write('Initial distance of elongation %s AA' % str(initial_dist) + '\n')
out.write('Step distance of pull-push intermediate steps %s AA' % str(Dl)+'\n')
out.close()

# check tinker commands (DEBUG for desired pulling speed)
opt = "dynamic "+'dummy'+" "+str(steps_thermalization)+" "+str(int_time)+" "\
    + str(dumps_time_therm)+" "+str(mode)+" "+str(temp_ener) \
    + " >> "+'outputxt'
os.system('echo \'This is the thermalization command\' >> %s ' % outputdat)
os.system('echo \'%s\' >> %s ' % (str(opt), outputdat))

copt = "dynamic "+'str(dummy)'+' '+str(steps)+" "+str(int_time)+" "+str(dump)\
    + " "+str(mode)+' '+str(temp_ener)+" >> "+'outputxt'

os.system('sleep 5s')
os.system('echo \'This is a generic dynamic intermediate step\' >> %s '
          % outputdat)
os.system('echo \'%s\' >> %s ' % (str(copt), outputdat))

os.system('echo \'Elongation speed requested is %s AA/ps \'>> %s '
          % (str(vel), outputdat))
actual = dl/(int(steps)*int_time*1e-3)
os.system('echo \'Actual pulling speed is %s AA/ps \' >> %s ' % (str(actual),
          outputdat))
os.system('echo \'Dumping %s structures per substep \' >> %s' % (str(ndumpstr),
          outputdat))
os.system('echo \'Dumping 1 structure every %.4g ps \' >> %s' % (dump,
          outputdat))
os.system('echo \'%.4g MD steps in one substep \' >> %s ' % (float(steps),
          outputdat))
os.system('echo \'Integration time for an MD step is %s fs \' >> %s' %
          (str(int_time), outputdat))

if abs(actual - vel) != 0:
    os.system('echo \'Error in pulling speed \' >> %s ' % outputdat)
    exit()
else:
    os.system('echo \'Checking speed: it is all good, keep simulating \' >> %s '
              % outputdat)

# define if tip and surface will be simulated or not
# tip = bool(False)
#  define the machine where it is running
# bluehive = bool(True)
# steno = bool(False)
#
# find carbon bound to endatoms
# if tip:
#     os.system('echo \'Tip is ACTIVE \' >> %s ' % outputdat)
#     # os.system('tinker2xyz_ %s' input_mol')
#     if bluehive:
#         os.system('obabel -itxyz %s -oxyz -O%s.xyz' % (input_mol, input_mol))
#         os.system('tinker2xyz %s' % input_mol)
#         os.system('tinker2xyz_ %s' % input_mol)
#     if steno:
#         os.system('tinker2xyz %s' % input_mol)
#         os.system('tinker2xyz_ %s' % input_mol)
#     moleculexyz = input_mol+'.xyz'
#     molecule = read('%s' % str(moleculexyz))
#     c1 = find_nearestC(molecule[int(argv[2])-1], molecule) + 1
#     c2 = find_nearestC(molecule[int(argv[3])-1], molecule) + 1
#     # fragments = Cluster(molecule)
#     fragment1 = fragments.find_connected(int(argv[3])-1)
#     fragment2 = fragments.find_connected(int(argv[2])-1)
#     list_frag1 = []
#     list_frag2 = []
#     for a in fragment1:
#         for b in molecule:
#             if a.position[0] == b.position[0] \
#              and a.position[1] == b.position[1] \
#              and a.position[2] == b.position[2]:
#                 list_frag1.append(b.index+1)
#     for a in fragment2:
#         for b in molecule:
#             if a.position[0] == b.position[0] \
#              and a.position[1] == b.position[1] \
#              and a.position[2] == b.position[2]:
#                 list_frag2.append(b.index+1)
# else:
#     os.system('echo \'Tip is NOT ACTIVE \' >> %s ' % outputdat)
# os.system('echo \'Tip is NOT ACTIVE \' >> %s ' % outputdat)

# Input molecule
input_molecule = open(argv[1], "r")

# create a key file with the same name of input xyz
name = str(argv[1])
dummy = name.replace(".xyz", "_dummy.xyz")
key = dummy.replace(".xyz", ".key")

# State the end atoms
endatom2 = str(argv[2])
endatom1 = str(argv[3])

# need to read the total number of atoms in the molecule
firstline = input_molecule.readline()
firstlinesp = firstline.split()
atomsnumb = int(firstlinesp[0])

# produce the dummy atoms number and the new total atom number
atoms = atomsnumb + 2
dummy1 = atomsnumb + 1
dummy2 = atoms
dummy3 = atoms + 1
dummy4 = atoms + 2
# if tip:
#     dummy5 = atoms + 3
#     dummy6 = atoms + 4

# now I need to read the coordinates of the endatoms
for line in input_molecule:
    linesp = line.split()
    if int(endatom1) == int(linesp[0]):
        dum1line = linesp
    elif int(endatom2) == int(linesp[0]):
        dum2line = linesp

input_molecule.close()
input_molecule = open(argv[1], "r")
# write the coordinates of dummy atoms
x1 = dum1line[2]
y1 = dum1line[3]
z1 = dum1line[4]
x2 = dum2line[2]
y2 = dum2line[3]
z2 = dum2line[4]
v1 = np.array([float(x1), float(y1), float(z1)])
v2 = np.array([float(x2), float(y2), float(z2)])
vd = v1-v2
# use this line to generate the initial condition for parallelization
v2 = (-(vd/mag(vd))*(initial_dist+0.00)) + v1
x2 = v2[0]
y2 = v2[1]
z2 = v2[2]
# translate the vector which will give the new coordinates to the dummy atoms
# the direction is chosen as the vector between the dummy atoms
v1t = (-(vd/mag(vd))*20.0) + v2
v2t = (-(vd/mag(vd))*22.0) + v2
# dummy atoms from the surface side
# if tip:
#     v3t = ((vd/mag(vd))*20.0) + v1
#     v4t = ((vd/mag(vd))*22.0) + v1

# create a new tinker file with dummy atoms and add the new
# lines of the dummy atons in a new tinker file
input_pull = open(str(dummy), "w")
# if tip:
#     input_pull.write(" "+"%s" % str(atoms+4)+"\n")
# else:
#     input_pull.write(" "+"%s" % str(atoms+2)+"\n")
input_pull.write(" "+"%s" % str(atoms+2)+"\n")

just1stline = input_molecule.readline()
for line in input_molecule:
    if line.strip():
        input_pull.write(line)

input_pull.write(" %s  Po   %s  %s  %s       200 " % (dummy1, x1, y1, z1,)+"\n")
input_pull.write(" %s  Po   %s  %s  %s       200 " % (dummy2, x2, y2, z2,)+"\n")
input_pull.write(" %s  Po   %s  %s  %s       200 " % (dummy3, v1t[0], v1t[1],
                                                      v1t[2])+"\n")
input_pull.write(" %s  Po   %s  %s  %s       200 " % (dummy4, v2t[0], v2t[1],
                                                      v2t[2])+"\n")
# if tip:
#     input_pull.write(" %s  Po   %s  %s  %s     200 " % (dummy5, v3t[0],
#                                                         v3t[1], v3t[2])+"\n")
#     input_pull.write(" %s  Po   %s  %s  %s     200 " % (dummy6, v4t[0],
#                                                         v4t[1], v4t[2])+"\n")

input_pull.close()

if thermalization:
    # now I create the key file making the dummy atoms inactive
    input_key = open(str(key), "w")
    input_key.write("""ARCHIVE
RANDOMSEED 123
THERMOSTAT NOSE-HOOVER
PARAMETERS /Users/alessap/tinker_hbonds/params/mm3_hbonds.prm
# PARAMETERS /Users/alessap/tinker_hbonds/params/mm3.prm
# PARAMETERS /kemi/alessap/bin/tinker/params/mm3_hbonds.prm
ATOM 200 200 Po \"Dummy Atom\" 0 0.0 0
angle     2 6   2      0.770       106.80                    #mm3
torsion 151 1   1   6  0.000 0.0 1 0.000 180.0 2 0.000 0.0 3 #mm3pro 6-1-1-8
torsion   1 1 151   3 -0.300 0.0 1 0.000 180.0 2 0.300 0.0 3 #mm3pro 1-1-8-2
torsion 151 1   1 151  0.000 0.0 1 0.000 180.0 2 0.000 0.0 3 #mm3pro 8-1-1-17
""")
    # if tip:
    #     input_key.write('INACTIVE ' + str(dummy1) + '  ' + str(dummy2) +
    #                     '  ' + str(dummy3) + '  ' + str(dummy4) + '  ' +
    #                     str(dummy5) + '   ' + str(dummy6) + "\n")
    # else:
    #     input_key.write('INACTIVE ' + str(dummy1) + '  ' + str(dummy2) +
    #                     '  ' + str(dummy3) + '  ' + str(dummy4) + "\n")
    input_key.write('INACTIVE ' + str(dummy1) + '  ' + str(dummy2) +
                    '  ' + str(dummy3) + '  ' + str(dummy4) + "\n")
    # this is the pulling extremity
    input_key.write('RESTRAIN-DISTANCE ' + str(endatom2) + ' ' + str(dummy2) +
                    ' ' + str(kcantilever) + ' 0.00 0.00' + '\n')
    input_key.write('RESTRAIN-ANGLE ' + str(dummy4) + ' ' + str(dummy3) +
                    ' ' + str(endatom2) + ' 10.00 180.0 180.0' + '\n')
#    if tip:
#        k_tip_1 = 0.5e-3 / len(list_frag1)
#        k_tip_2 = 1e-3 / len(list_frag1)
#        k_tip_3 = 2e-3 / len(list_frag1)
#        for ind in list_frag2:
#            input_key.write('RESTRAIN-ANGLE ' + str(dummy4) + ' ' +
#                            str(dummy2) + ' ' + str(ind) + ' ' + str(k_tip_1) +
#                            ' 90.0 270.0 ' + '\n')
#            input_key.write('RESTRAIN-ANGLE ' + str(dummy4) + ' ' +
#                            str(dummy2) + ' ' + str(ind) + ' ' +
#                            str(k_tip_2)+' 80.0 280.0 '+'\n')
#            input_key.write('RESTRAIN-ANGLE ' + str(dummy4) + ' ' +
#                            str(dummy2) + ' ' + str(ind) + ' ' +
#                            str(k_tip_3) + ' 70.0 290.0 ' + '\n')
#        for ind in list_frag1:
#            input_key.write('RESTRAIN-ANGLE '+str(dummy4)+' '+str(dummy2) +
#                            ' '+str(ind)+' '+str(k_tip_1)+' 90.0 270.0 '+'\n')
#            input_key.write('RESTRAIN-ANGLE '+str(dummy4)+' '+str(dummy2)+' ' +
#                            str(ind)+' '+str(k_tip_2)+' 80.0 280.0 '+'\n')
#            input_key.write('RESTRAIN-ANGLE '+str(dummy4)+' '+str(dummy2)+' ' +
#                            str(ind)+' '+str(k_tip_3)+' 70.0 290.0 '+'\n')

    # this is the atom stuck onto the surface
    input_key.write('RESTRAIN-DISTANCE '+str(endatom1)+' '+str(dummy1) +
                    ' 2000.00 0.00 0.00'+'\n')
#    if tip:
#        k_surface_1 = 0.5e-3 / len(list_frag1)
#        k_surface_2 = 1e-3 / len(list_frag1)
#        k_surface_3 = 2e-3 / len(list_frag1)
#        for ind in list_frag1:
#            input_key.write('RESTRAIN-ANGLE '+str(dummy6)+' '+str(dummy1)+' ' +
#                            str(ind) + ' ' + str(k_surface_1) +
#                            ' 90.0 270.0 ' + '\n')
#            input_key.write('RESTRAIN-ANGLE '+str(dummy6)+' '+str(dummy1)+' ' +
#                            str(ind) + ' '+str(k_surface_2) +
#                            ' 80.0 280.0 '+'\n')
#            input_key.write('RESTRAIN-ANGLE '+str(dummy6)+' '+str(dummy1)+' ' +
#                            str(ind)+' '+str(k_surface_3) +
#                            ' 70.0 290.0 '+'\n')
#        for ind in list_frag2:
#            input_key.write('RESTRAIN-ANGLE '+str(dummy6)+' '+str(dummy1)+' ' +
#                            str(ind)+' '+str(k_surface_1) +
#                            ' 90.0 270.0 '+'\n')
#            input_key.write('RESTRAIN-ANGLE '+str(dummy6)+' '+str(dummy1)+' ' +
#                            str(ind)+' '+str(k_surface_2) +
#                            ' 80.0 280.0 '+'\n')
#            input_key.write('RESTRAIN-ANGLE '+str(dummy6)+' '+str(dummy1)+' ' +
#                            str(ind)+' '+str(k_surface_3) +
#                            ' 70.0 290.0 '+'\n')

    input_key.close()
    input_molecule.close()

    # optimize the geometry with the dummy atoms
    # outputxt="%s_output.txt" % name.replace(".xyz","")
    outputxt = "md_output.txt"
    output_opt = open(str(outputxt), "w")
    theline = "Optimization of the system with dummy atoms constrainsts  - %s"\
        % name
    output_opt.write(theline+"\n")
    opt = "minimize "+dummy+" 1e-1  "+" >> "+outputxt
    output_opt.write(opt+"\n")
    os.system(opt)
    os.system('mv '+dummy+'_2 '+dummy)

    opt = "newton "+dummy+" A A 1e-10  "+" >> "+outputxt
    output_opt.write(opt+"\n")
    os.system(opt)
    output_opt.close()
    os.system('mv '+dummy+'_2 '+dummy)

    # thermalization of the initial geometry
    # outputxt="%s_output.txt" % name.replace(".xyz","")
    output_opt = open(str(outputxt), "w")
    theline = "Thermalisation of the molecular system with end atoms \
           constraints - %s" % name
    output_opt.write(theline+"\n")
    # if tip:
    #     for incrt in range(10, 300, 10):
    #         opt = "dynamic "+dummy+" 30000 1.0 30 2 "+str(incrt)+" >> "+outputxt
    #         os.system(opt)
    #         os.system('sleep 5s')
    #     os.system('cp '+dummy+' partly_thermalized_t.xyz')
    opt = "dynamic " + dummy + " " + str(steps_thermalization) + " " + \
          str(int_time) + " " + str(dumps_time_therm) + " " + str(mode) + \
          " " + str(temp_ener)+" >> "+outputxt
    os.system('sleep 5s')
    output_opt.write(opt+"\n")
    os.system('echo \'Thermalizing for %s  ns \' >> %s ' % (str(vel),
              outputdat))
    os.system('echo \'Thermalizion started at ' +
              str(datetime.datetime.now().isoformat())+'  \'  >> %s '
              % outputdat)
    # os.system('echo Sun Sep 25 11:34:22 CEST 2016 >> %s ' % outputdat)
    os.system(opt)
    output_opt.close()
    os.system('echo \'Thermalizion finished at ' +
              str(datetime.datetime.now().isoformat())+'  \' >> %s '
              % outputdat)
    # os.system('echo Sun Sep 25 11:34:22 CEST 2016 >> %s ' % outputdat)
    os.system('cp '+dummy+' thermalized_t.xyz')
    os.system('cp ' + str(dummy).replace(".xyz", ".dyn") +
              '  thermalized_t.dyn')
# else:
#     outputxt = "md_output.txt"
#     output_opt = open(str(outputxt), "w")
#     theline = "Optimization of the system with dummy atoms constrainsts  \
#               - %s" % name
#     output_opt.write(theline+"\n")
#     opt = "minimize "+dummy+" 1e-1  "+" >> "+outputxt
#     output_opt.write(opt+"\n")
#     # os.system(opt)
#     # os.system('mv '+dummy+'_2 '+dummy)
#     opt = "newton "+dummy+" A A 1e-10  "+" >> "+outputxt
#     output_opt.write(opt+"\n")
#     # os.system(opt)
#     output_opt.close()
#     # os.system('mv '+dummy+'_2 '+dummy)
#     # thermalization of the initial geometry
#     # outputxt="%s_output.txt" % name.replace(".xyz","")
#     output_opt = open(str(outputxt), "w")
#     theline = "Thermalisation of the molecular system with end atoms \
#               constraints - %s" % name
#     output_opt.write(theline+"\n")
#     opt = "dynamic "+dummy+" "+str(steps_thermalization)+" "+str(int_time) + \
#           " "+str(dumps_time_therm)+" "+str(mode)+" "+str(temp_ener) + \
#           " >> "+outputxt
#     os.system('sleep 5s')
#     output_opt.write(opt+"\n")
#     os.system('echo \'Thermalizing for 5.0  ns \' >> %s ' % outputdat)
#     os.system('echo \'Thermalizion started at ' +
#               str(datetime.datetime.now().isoformat()) +
#               '  \'  >> %s ' % outputdat)
#     # os.system('echo Sun Sep 25 11:34:22 CEST 2016 >> %s ' % outputdat)
#     # os.system(opt)
#     output_opt.close()
#     os.system('echo \'Thermalizion finished at ' +
#               str(datetime.datetime.now().isoformat()) +
#               '  \' >> %s ' % outputdat)
#     # os.system('echo Sun Sep 25 11:34:22 CEST 2016 >> %s ' % outputdat)
#     # os.system('cp '+dummy+' thermalized_t.xyz')
#     # os.system('cp '+str(dummy).replace(".xyz",".dyn")+'  thermalized_t.dyn')
#     #########
#     lastdir = '/kemi/alessap/exciton/pulling/PDI-PDI/CH2_2/..\
#               /step99_9.8to9.9/'
#     os.system('cp '+str(lastdir)+'thermalized_t.xyz.gz '+str(lastdir) +
#               'thermalized_t.dyn.gz '+str(lastdir) +
#               str(dummy).replace(".xyz", ".key")+'.gz . ')
#     os.system('gzip -d thermalized_t.xyz.gz thermalized_t.dyn.gz  ' +
#               str(dummy).replace(".xyz", ".key")+'.gz   ')
#     os.system('cp '+str(lastdir)+'thermalized_t.xyz.gz '+str(lastdir) +
#               'thermalized_t.dyn.gz '+str(lastdir) +
#               str(dummy).replace(".xyz", ".key")+'.gz . ')
#     os.system('cp thermalized_t.xyz '+str(dummy))
#     os.system('cp thermalized_t.xyz '+str(dummy).replace(".xyz", ".arc"))
#     os.system('cp thermalized_t.dyn '+str(dummy).replace(".xyz", ".dyn"))

# removing the first snapshot
# os.system('tail '+'-'+str(int(dummy6)+1)+' thermalized_t.xyz >'+dummy)

e = "D+00"

# print "This was a test "
# exit()

##############################################################################
#                                                                            #
#                 Ready to do the elongation-contraction cycle               #
#                                                                            #
##############################################################################
dltot = dl
t = 0
structures = 0
force = open('force.dat', 'a')
os.system('echo \'Starting ELONGATION \' >> %s ' % outputdat)

for x in range(n):
    structures = structures+1
    # read the coordinates of the dummy atoms from the dyn file
    input = open(str(dummy).replace(".xyz", ".dyn"), "r")
    lines = input.readlines()
    input.close()
    dum1line = lines[dummy1+5]  # 5 is the number of line at the top of the dyn
    dum2line = lines[dummy2+5]
    dum3line = lines[dummy3+5]
    dum4line = lines[dummy4+5]
    # if tip:
    #     dum5line = lines[dummy5+5]
    #     dum6line = lines[dummy6+5]
    v1 = np.array([float(x.replace("D", "e")) for x in dum1line.split()])
    v2 = np.array([float(x.replace("D", "e")) for x in dum2line.split()])
    v3 = np.array([float(x.replace("D", "e")) for x in dum3line.split()])
    v4 = np.array([float(x.replace("D", "e")) for x in dum4line.split()])
    # if tip:
    #     v5 = np.array([float(x.replace("D", "e")) for x in dum5line.split()])
    #     v6 = np.array([float(x.replace("D", "e")) for x in dum6line.split()])
    # pull move the dummy atom(s)
    vd = v1-v2
    L = mag(vd)
    # translate the vector which will give the new coord. to the dummy atoms
    # the direction is chosen as the vector between the dummy atoms
    v2t = (-(vd/mag(vd))*dl)+v2
    v3t = (-(vd/mag(vd))*dl)+v3
    v4t = (-(vd/mag(vd))*dl)+v4
    v1t = v1    # atom stuck to the surface
    # if tip:
    #     v5t = v5    # atom stuck to the surface
    #     v6t = v6    # atom stuck to the surface
    # now I rewrite the input file for the elongation
    dum1linet = "    %s%s     %s%s    %s%s  " % (format(v1t[0], '.16f'), e,
                                                 format(v1t[1], '.16f'), e,
                                                 format(v1t[2], '.16f'), e,)
    dum2linet = "    %s%s     %s%s    %s%s  " % (format(v2t[0], '.16f'), e,
                                                 format(v2t[1], '.16f'), e,
                                                 format(v2t[2], '.16f'), e,)
    dum3linet = "    %s%s     %s%s    %s%s  " % (format(v3t[0], '.16f'), e,
                                                 format(v3t[1], '.16f'), e,
                                                 format(v3t[2], '.16f'), e,)
    dum4linet = "    %s%s     %s%s    %s%s  " % (format(v4t[0], '.16f'), e,
                                                 format(v4t[1], '.16f'), e,
                                                 format(v4t[2], '.16f'), e,)
    # if tip:
    #     dum5linet = "  %s%s     %s%s    %s%s  " % (format(v5t[0], '.16f'), e,
    #                                                format(v5t[1], '.16f'), e,
    #                                                format(v5t[2], '.16f'), e,)
    #     dum6linet = "  %s%s     %s%s    %s%s  " % (format(v6t[0], '.16f'), e,
    #                                                format(v6t[1], '.16f'), e,
    #                                                format(v6t[2], '.16f'), e,)
    # read the number of lines in dyn file
    with open(str(dummy).replace(".xyz", ".dyn"), 'r') as f:
        flength = len(f.readlines())
    input = open(str(dummy).replace(".xyz", ".dyn"), "w")
    for line in range(dummy1+5):
        input.write(lines[line])
    # if tip:
    #     input.write(dum1linet+"\n"+dum2linet+"\n"+dum3linet+"\n"+dum4linet +
    #                 "\n" + dum5linet+"\n"+dum6linet+"\n")
    #     for line in range(dummy6+5+1, flength):
    #         input.write(lines[line])
    # else:
    #   input.write(dum1linet+"\n"+dum2linet+"\n"+dum3linet+"\n"+dum4linet+"\n")
    #     for line in range(dummy4+5+1, flength):
    #         input.write(lines[line])
    input.write(dum1linet+"\n"+dum2linet+"\n"+dum3linet+"\n"+dum4linet+"\n")
    for line in range(dummy4+5+1, flength):
        input.write(lines[line])
    input.close()
    # run the  MD simulation
    copt = "dynamic " + str(dummy) + '  '+str(steps) + " " + str(int_time) + \
           " " + str(dump) + " " + str(mode) + ' ' + str(temp_ener) + \
           " >> " + outputxt
    os.system(copt)
    # now I need to calculate all the forces during the MD step run before
    # ndumpstr : structure dumped per md step
    # sed the last ndumped structure in a temporary arc file
    # if tip:
    #     linees = int((int(t)*int(ndumpstr)*int(dummy6 + 1)) + int(dummy6 + 1))
    #     os.system("sed '1,"+str(linees)+"d' " +
    #               str(str(dummy).replace('xyz', 'arc'))+" > temp.arc")
    # else:
    #     linees = int((int(t)*int(ndumpstr)*int(dummy4+1))+int(dummy4+1))
    #     os.system("sed '1,"+str(linees)+"d' " +
    #               str(str(dummy).replace('xyz', 'arc'))+" > temp.arc")
    linees = int((int(t)*int(ndumpstr)*int(dummy4+1))+int(dummy4+1))
    os.system("sed '1,"+str(linees)+"d' " +
              str(str(dummy).replace('xyz', 'arc'))+" > temp.arc")
    # loop over ndumped structure
    temp = open('temp.arc', 'r')
    templ = temp.readlines()
    for i in range(int(ndumpstr)):
        # if tip:
        #     line1 = templ[((int(i)*int(dummy6+1))+int(endatom1))]
        #     line2 = templ[((int(i)*int(dummy6+1))+int(endatom2))]
        # else:
        #     line1 = templ[((int(i)*int(dummy4+1))+int(endatom1))]
        #     line2 = templ[((int(i)*int(dummy4+1))+int(endatom2))]
        line1 = templ[((int(i)*int(dummy4+1))+int(endatom1))]
        line2 = templ[((int(i)*int(dummy4+1))+int(endatom2))]
        linez1 = line1.split()
        linez2 = line2.split()
        z1 = np.array([float(linez1[2]), float(linez1[3]), float(linez1[4])])
        z2 = np.array([float(linez2[2]), float(linez2[3]), float(linez2[4])])
        zeta = z1-z2
        Z = mag(zeta)
        # cantilever force constant (in kcal/A^2)
        # (1kcal/A^2mol = 69.5 pN/A=0.695N/m)
        F = -kcantilever*(float(Z)-float(L))*69.5
        # force.write('%.0f %.4f %.8f %.8f %.8f %.8f' %
        # (n,(dltot+initial_dist),Z,L,F,Z-L)+'\n')
        force.write('0   1   %.8f %.8f %.8f ' % (Z, L, F)+'\n')
    temp.close()
    t = t+1  # counter of cycle
    # os.system('echo \'Step %.0f/%.0f completed\' '
    # +str(datetime.datetime.now().isoformat())+'  >> %s ' % (t,n*2,outputdat))
    # os.system('echo \'Thermalizion finished at '
    # + str(datetime.datetime.now().isoformat())+'  \' >> %s ' %  outputdat)
    # os.system('echo \'Step %.0f/%.0f completed '
    # +str(datetime.datetime.now().isoformat())
    # +'  \' >> %s ' % (t,n*2,outputdat))
    # os.system('echo Sun Sep 25 11:34:22 CEST 2016 >> %s ' % outputdat)
    os.system('echo Step %s / %s completed  at %s >> %s ' % (str(t),
              str(n*2), str(datetime.datetime.now().isoformat()), outputdat))
    dltot = dl + dltot

tt = 0
dltot = dl
Dl1 = Dl+initial_dist
os.system('echo \'Starting CONTRACTION \' >> %s ' % outputdat)
for x in range(n):
    structures = structures+1
    # read the coordinates of the dummy atoms from the dyn file
    input = open(str(dummy).replace(".xyz", ".dyn"), "r")
    lines = input.readlines()
    input.close()
    dum1line = lines[dummy1+5]  # 5 is the number of line at the top of the dyn
    dum2line = lines[dummy2+5]
    dum3line = lines[dummy3+5]
    dum4line = lines[dummy4+5]
    # if tip:
    #     dum5line = lines[dummy5+5]
    #     dum6line = lines[dummy6+5]
    v1 = np.array([float(x.replace("D", "e")) for x in dum1line.split()])
    v2 = np.array([float(x.replace("D", "e")) for x in dum2line.split()])
    v3 = np.array([float(x.replace("D", "e")) for x in dum3line.split()])
    v4 = np.array([float(x.replace("D", "e")) for x in dum4line.split()])
    # if tip:
    #     v5 = np.array([float(x.replace("D", "e")) for x in dum5line.split()])
    #     v6 = np.array([float(x.replace("D", "e")) for x in dum6line.split()])
    # pull move the dummy atom(s)
    vd = v1-v2
    L = mag(vd)
    # translate the vector which will give the new coord. to the dummy atoms
    # the direction is chosen as the vector between the dummy atoms
    v2t = (+(vd/mag(vd))*dl)+v2
    v3t = (+(vd/mag(vd))*dl)+v3
    v4t = (+(vd/mag(vd))*dl)+v4
    v1t = v1    # atom stuck to the surface
    # if tip:
    #     v5t = v5    # atom stuck to the surface
    #     v6t = v6    # atom stuck to the surface
    # now I rewrite the input file for the elongation
    dum1linet = "    %s%s     %s%s    %s%s  " % (format(v1t[0], '.16f'), e,
                                                 format(v1t[1], '.16f'), e,
                                                 format(v1t[2], '.16f'), e,)
    dum2linet = "    %s%s     %s%s    %s%s  " % (format(v2t[0], '.16f'), e,
                                                 format(v2t[1], '.16f'), e,
                                                 format(v2t[2], '.16f'), e,)
    dum3linet = "    %s%s     %s%s    %s%s  " % (format(v3t[0], '.16f'), e,
                                                 format(v3t[1], '.16f'), e,
                                                 format(v3t[2], '.16f'), e,)
    dum4linet = "    %s%s     %s%s    %s%s  " % (format(v4t[0], '.16f'), e,
                                                 format(v4t[1], '.16f'), e,
                                                 format(v4t[2], '.16f'), e,)
    # if tip:
    #     dum5linet = "    %s%s     %s%s    %s%s  " % (format(v5t[0], '.16f'), e,
    #                                                  format(v5t[1], '.16f'), e,
    #                                                  format(v5t[2], '.16f'), e,)
    #     dum6linet = "    %s%s     %s%s    %s%s  " % (format(v6t[0], '.16f'), e,
    #                                                  format(v6t[1], '.16f'), e,
    #                                                  format(v6t[2], '.16f'), e,)
    # read the number of lines in dyn file
    with open(str(dummy).replace(".xyz", ".dyn"), 'r') as f:
        flength = len(f.readlines())
    input = open(str(dummy).replace(".xyz", ".dyn"), "w")
    for line in range(dummy1+5):
        input.write(lines[line])
    # if tip:
    #     input.write(dum1linet+"\n"+dum2linet+"\n"+dum3linet+"\n"+dum4linet +
    #                 "\n" + dum5linet+"\n"+dum6linet+"\n")
    #     for line in range(dummy6+5+1, flength):
    #         input.write(lines[line])
    # else:
    #     input.write(dum1linet+"\n"+dum2linet+"\n"+dum3linet+"\n"+dum4linet+"\n")
    #     for line in range(dummy4+5+1, flength):
    #         input.write(lines[line])
    input.write(dum1linet+"\n"+dum2linet+"\n"+dum3linet+"\n"+dum4linet+"\n")
    for line in range(dummy4+5+1, flength):
        input.write(lines[line])
    input.close()
    # run the  MD simulation
    copt = "dynamic " + str(dummy) + '  '+str(steps) + " " + str(int_time) + \
           " " + str(dump) + " " + str(mode) + ' ' + str(temp_ener) + \
           " >> " + outputxt
    os.system(copt)
    # now I need to calculate all the forces during the MD step run before
    # ndumpstr : structure dumped per md step
    # sed the last ndumped structure in a temporary arc file
    # if tip:
    #     linees = int((int(t)*int(ndumpstr)*int(dummy6 + 1)) + int(dummy6 + 1))
    #     os.system("sed '1,"+str(linees)+"d' " +
    #               str(str(dummy).replace('xyz', 'arc'))+" > temp.arc")
    # else:
    #     linees = int((int(t)*int(ndumpstr)*int(dummy4+1))+int(dummy4+1))
    #     os.system("sed '1,"+str(linees)+"d' " +
    #               str(str(dummy).replace('xyz', 'arc'))+" > temp.arc")
    linees = int((int(t)*int(ndumpstr)*int(dummy4+1))+int(dummy4+1))
    os.system("sed '1,"+str(linees)+"d' " +
              str(str(dummy).replace('xyz', 'arc'))+" > temp.arc")
    # loop over ndumped structure
    temp = open('temp.arc', 'r')
    templ = temp.readlines()
    for i in range(int(ndumpstr)):
        # if tip:
        #     line1 = templ[((int(i)*int(dummy6+1))+int(endatom1))]
        #     line2 = templ[((int(i)*int(dummy6+1))+int(endatom2))]
        # else:
        #     line1 = templ[((int(i)*int(dummy4+1))+int(endatom1))]
        #     line2 = templ[((int(i)*int(dummy4+1))+int(endatom2))]
        line1 = templ[((int(i)*int(dummy4+1))+int(endatom1))]
        line2 = templ[((int(i)*int(dummy4+1))+int(endatom2))]
        linez1 = line1.split()
        linez2 = line2.split()
        z1 = np.array([float(linez1[2]), float(linez1[3]), float(linez1[4])])
        z2 = np.array([float(linez2[2]), float(linez2[3]), float(linez2[4])])
        zeta = z1-z2
        Z = mag(zeta)
        # cantilever force constant (in kcal/A^2)
        # (1kcal/A^2mol = 69.5 pN/A=0.695N/m)
        F = -kcantilever*(float(Z)-float(L))*69.5
        # force.write('%.0f %.4f %.8f %.8f %.8f %.8f' %
        # (n,(dltot+initial_dist),Z,L,F,Z-L)+'\n')
        force.write('0   1   %.8f %.8f %.8f ' % (Z, L, F)+'\n')
    temp.close()
#    structures=structures+1
#    # Read the coordinates of the dummy atoms from the dyn file
#    input=open(str(dummy).replace(".xyz",".dyn"), "r")
#    lines=input.readlines()
#    input.close()
#    dum1line=lines[dummy1+5] # 5 is the number of line at the top of the dyn
#    dum2line=lines[dummy2+5]
#    dum3line=lines[dummy3+5]
#    dum4line=lines[dummy4+5]
#    if tip:
#      dum5line=lines[dummy5+5]
#      dum6line=lines[dummy6+5]
#    v1 =  np.array([float(x.replace("D","e")) for x in dum1line.split()])
#    v2 =  np.array([float(x.replace("D","e")) for x in dum2line.split()])
#    v3 =  np.array([float(x.replace("D","e")) for x in dum3line.split()])
#    v4 =  np.array([float(x.replace("D","e")) for x in dum4line.split()])
#    if tip:
#      v5 =  np.array([float(x.replace("D","e")) for x in dum5line.split()])
#      v6 =  np.array([float(x.replace("D","e")) for x in dum6line.split()])
#    #b move the dummy atoms
#    vd = v1-v2
#    L = mag(vd)
#    # translate the vector which will give the new coords to the dummy atoms
#    # the direction is chosen as the vector between the dummy atoms
#    v2t=((vd/mag(vd))*dl)+v2
#    v3t=((vd/mag(vd))*dl)+v3
#    v4t=((vd/mag(vd))*dl)+v4
#    v1t=v1 # atom stuck to the surface
#    if tip:
#      v5t=v5    # atom stuck to the surface
#      v6t=v6    # atom stuck to the surface
#    # now I rewrite the input file for the elongation
#    dum1linet= "    %s%s     %s%s    %s%s  " % (format(v1t[0],'.16f'),e,
# format(v1t[1],'.16f'),e,format(v1t[2],'.16f'),e,)
#    dum2linet= "    %s%s     %s%s    %s%s  " % (format(v2t[0],'.16f'),e,
# format(v2t[1],'.16f'),e,format(v2t[2],'.16f'),e,)
#    dum3linet= "    %s%s     %s%s    %s%s  " % (format(v3t[0],'.16f'),e,
# format(v3t[1],'.16f'),e,format(v3t[2],'.16f'),e,)
#    dum4linet= "    %s%s     %s%s    %s%s  " % (format(v4t[0],'.16f'),e,
# format(v4t[1],'.16f'),e,format(v4t[2],'.16f'),e,)
#    if tip:
#      dum5linet= "    %s%s     %s%s    %s%s  " % (format(v5t[0],'.16f'),e,
# format(v5t[1],'.16f'),e,format(v5t[2],'.16f'),e,)
#      dum6linet= "    %s%s     %s%s    %s%s  " % (format(v6t[0],'.16f'),e,
# format(v6t[1],'.16f'),e,format(v6t[2],'.16f'),e,)
#    #read the number of lines in dyn file
#    with open(str(dummy).replace(".xyz",".dyn"), 'r') as f:
#      flength = len(f.readlines())
#    input=open(str(dummy).replace(".xyz",".dyn"), "w")
#    for line in range(dummy1+5):
#      input.write(lines[line])
#    if tip:
#      input.write(dum1linet+"\n"+dum2linet+"\n"+dum3linet+"\n"+dum4linet+"\n"
# +dum5linet+"\n"+dum6linet+"\n")
#      for line in range(dummy6+5+1, flength):
#        input.write(lines[line])
#    else:
#      input.write(dum1linet+"\n"+dum2linet+"\n"+dum3linet+"\n"+dum4linet+"\n")
#      for line in range(dummy4+5+1, flength):
#        input.write(lines[line])
#    input.close()
#    #
#    #run the  MD simulation
#    copt = "dynamic "+str(dummy)+'  '+str(steps)+" "+str(int_time)+" "+
# str(dump)+" "+str(mode)+' '+str(temp_ener)+" >> "+outputxt
#    os.system(copt)
#    # now I need to calculate all the forces during the MD step run before
#    # ndumpstr : structure dumped per md step
#    # sed the last ndumped structure in a temporary arc file
#    if tip:
#      linees=int((int(t+tt)*int(ndumpstr)*int(dummy6+1))+int(dummy6+1))
#      os.system("sed '1,"+str(linees)+"d' "+str(str(dummy).replace('xyz',#
# 'arc'))+" > temp.arc")
#    else:
#      linees=int((int(t+tt)*int(ndumpstr)*int(dummy4+1))+int(dummy4+1))
#      os.system("sed '1,"+str(linees)+"d' "+str(str(dummy).replace('xyz',
# 'arc'))+" > temp.arc")
#      #os.system("sed '1,"+str(dummy4+1)+"d' "+str(str(dummy).replace('xyz',
# 'arc'))+" > temp.arc")
#    # loop over ndumped structure
#    temp=open('temp.arc','r')
#    templ=temp.readlines()
#    for i in range(int(ndumpstr)):
#      if tip:
#        line1=templ[((int(i)*int(dummy6+1))+int(endatom1))]
#        line2=templ[((int(i)*int(dummy6+1))+int(endatom2))]
#      else:
#        line1=templ[((int(i)*int(dummy4+1))+int(endatom1))]
#        line2=templ[((int(i)*int(dummy4+1))+int(endatom2))]
#      linez1=line1.split()
#      linez2=line2.split()
#      z1 =  np.array([float(linez1[2]),float(linez1[3]),float(linez1[4])])
#      z2 =  np.array([float(linez2[2]),float(linez2[3]),float(linez2[4])])
#      zeta=z1-z2
#      Z = mag(zeta)
#      # cantilever force constant (in kcal/A^2)
# (1kcal/A^2mol = 69.5 pN/A=0.695N/m)
#      F=-kcantilever*(float(Z)-float(L))*69.5
#      force.write('1   1   %.8f %.8f %.8f ' % (Z,L,F)+'\n')
#    temp.close()
    tt = tt+1  # counter of cycle
    # os.system('echo \'Step %.0f/%.0f completed ' % (tt+t, n*2) +
    #           str(datetime.datetime.now().isoformat()) +
    #           '  \' >> %s ' % (outputdat))
    # os.system('echo \'Step %.0f/%.0f completed\' ' % (tt+t, n*2) +
    #           str(datetime.datetime.now().isoformat()) +
    #           '  >> %s ' % (outputdat))
    # os.system('echo \'Step %.0f/%.0f completed\'   >> %s ' % (tt+t, n*2,
    #           outputdat))
    #          (tt+t, n*2, outputdat))
    # os.system('echo Sun Sep 25 11:34:22 CEST 2016 >> %s ' % outputdat)
    os.system('echo Step %s / %s completed  at %s >> %s ' %
              (str(tt+t), str(n*2),
               str(datetime.datetime.now().isoformat()), outputdat))
    dltot = dl + dltot

##############################################################################
#                                                                            #
#                        END of the pulling cycle                            #
#                                                                            #
##############################################################################
force.close()
os.system('echo \'DONE\' >> %s ' % outputdat)
exit()
