#!/auto/vestec1-elixir/home/manishkumar/.conda/envs/kpython310/bin/python3.1
import sys

fo= open("cas_orbital.molden.input",'w')
line_num = 0
list_of_res = []
basis_data = []


def skipline(n):
    for j in range(n):
        readobj.readline()

with open(sys.argv[1],'r') as readobj:
    for line in readobj:
        line_num+=1
        if "Sym" in line:
            list_of_res.append(line_num)
        if "[MO]" in line:
            basis_data.append(line_num)

last_row = list_of_res[int(sys.argv[3])+1]-list_of_res[int(sys.argv[2])]


# print(list_of_res[8])
# print(basis_data)
with open(sys.argv[1],'r') as readobj:
    x = readobj.tell()
    for i in range(int(basis_data[0])):
        fo.write(readobj.readline())
    readobj.seek(x)
    skipline(list_of_res[int(sys.argv[2])]-1)
    for j in range(int(last_row)):
        fo.write(readobj.readline())

