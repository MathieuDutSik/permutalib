import os;
import sys;

if len(sys.argv) != 3:
    print("convert [FileI] [FileO]")
    os.sys.exit(1)


FileI = sys.argv[1];
FileO = sys.argv[2];

ifile = open(FileI, 'r')
l_lines = []
for e_line in ifile:
    l_lines.append(e_line)


f = open(FileO, 'w')
n_group = int(l_lines[0])
pos = 1
f.write(str(n_group) + '\n')
for i_group in range(n_group):
    desc_line = l_lines[pos]
    pos += 1
    U = desc_line.split(" ")
    n_gen = int(U[0])
    n_act = int(U[1])
    f.write(str(n_act) + " " + str(n_gen) + "\n")
    for i in range(n_gen):
        gen_line = l_lines[pos]
        pos += 1
        f.write(gen_line)

f.close()
