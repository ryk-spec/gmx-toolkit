# ambdat2gmx.py M.Yoneya 30.06.2015
import sys


cal = 4.184
status = "initial"

# 変数の定義
mass = {}

status = "initial"

with open(sys.argv[1], "r") as infile:        
    for line in infile:
        line = line.strip()
        if line == "":
            continue
        if status == "initial" and "MASS" in line:
            status = "MASS" 
        elif status == "MASS" and line != "BOND":
            mass[line.split()[0]] = float(line.split()[1])
        elif status == "MASS" and line == "BOND":
            status = "BOND"
        elif status == "BOND" and line != "ANGLE":
            parts=line.replace("-"," ").split()
            print(
                f"{parts[0]:<10}{parts[1]:<10}{1:2}{0.1*float(parts[3]):>10.5f}{100*2*cal*float(parts[2]):>12.1f}"
            )
        elif status == "BOND" and line == "ANGLE":
            status = "ANGLE"
        elif status == "ANGLE" and line != "DIHE":
            parts = line.replace("-", " ").split()
            func = 1
            line = f"{parts[0]:<10}{parts[1]:<10}{parts[2]:<10}{func:2}{float(parts[4]):>10.3f}{2*cal*float(parts[3]):>10.3f}"
            print(
                line
            )
        elif status == "ANGLE" and line == "DIHE":
            status = "DIHE"
        elif status == "DIHE" and line != "NONBON":
            parts = line.replace("-", " ").split()
            func = 9
            line = f"{parts[0]:<10}{parts[1]:<10}{parts[2]:<10}{parts[3]:<10}{func:2}{float(parts[6]):>8.1f}{cal*float(parts[5])/float(parts[4]):>12.5f}{abs(float(parts[7])):>5.0f} "
            print(line)
        elif status == "DIHE" and line == "NONBON":
            status = "nonbonded"
        elif status == "nonbonded" and line != "NONBON":
            parts = line.replace("-", " ").split()
            if parts[0][:2] == "Na":
                atnum = 11
            elif parts[0][:2] in ["cl", "Cl"]:
                atnum = 17
            elif parts[0][:2] == "C0":
                atnum = 20
            elif parts[0][:2] in ["br", "Br"]:
                atnum = 35
            elif parts[0][:2] == "Cs":
                atnum = 55
            elif parts[0][:2] == "IM":
                atnum = 17
            elif parts[0][:2] == "IB":
                atnum = 0
            elif parts[0][:1] in ["h", "H"]:
                atnum = 1
            elif parts[0][:1] in ["c", "C"]:
                atnum = 6
            elif parts[0][:1] in ["n", "N"]:
                atnum = 7
            elif parts[0][:1] in ["o", "O"]:
                atnum = 8
            elif parts[0][:1] in ["f", "F"]:
                atnum = 9
            elif parts[0][:1] in ["p", "P"]:
                atnum = 15
            elif parts[0][:1] in ["s", "S"]:
                atnum = 16
            elif parts[0][:1] in ["i", "I"]:
                atnum = 53
            else:
                atnum = 0
            line = f"{parts[0]:<10}{atnum:>3}{mass[parts[0]]:>8.3f}  0.0  A {float(parts[1])*0.1*(2**(5/6)):.5e} {cal*float(parts[2]):.5e}"
            print(line)

