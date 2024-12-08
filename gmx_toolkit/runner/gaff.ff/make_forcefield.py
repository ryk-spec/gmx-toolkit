# ambdat2gmx.py M.Yoneya 30.06.2015
import sys


cal = 4.184
status = "initial"

# 変数の定義
mass = {}

# ファイル操作
with open("forcefield.doc", "w") as docfile, \
     open("atomtypes.atp", "w") as atpfile, \
     open("ffnonbonded.itp", "w") as nonfile, \
     open("ffbonded.itp", "w") as bonfile1, \
     open("ffbonded2.itp", "w") as bonfile2, \
     open("forcefield.itp", "w") as fffile:
    
    # ファイルの初期設定
    # docfile.write("#define _FF_AMBER\n")
    # docfile.write("#define _FF_AMBERGENERAL\n")
    # docfile.write("\n")
    fffile.write("[ defaults ]\n")
    fffile.write("; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n")
    fffile.write("1               2               yes             0.5     0.8333\n")
    fffile.write("\n")
    fffile.write("#include \"ffnonbonded.itp\"\n")
    fffile.write("#include \"ffbonded.itp\"\n")
    fffile.write("#include \"ffbonded2.itp\"\n")
    
    # ファイルのオープン
    status = "initial"
    with open(sys.argv[1], "r") as infile:
        for line in infile:
            line = line.strip()
            if line == "":
                if status == "masstypes":
                    status = "hydrophilic"
                    bonfile1.write("[ bondtypes ]\n")
                    bonfile1.write("; i    j  func       b0          kb\n")
                elif status == "bondtypes":
                    status = "angletypes"
                    bonfile1.write(f"{line}\n")
                    bonfile1.write("[ angletypes ]\n")
                    bonfile1.write(";  i    j    k  func       th0       cth\n")
                elif status == "angletypes":
                    status = "propertypes"
                    bonfile1.write(f"{line}\n")
                    bonfile2.write("[ dihedraltypes ] ; proper stuff\n")
                    bonfile2.write(";i  j   k  l	 func      phase      kd      pn\n")
                elif status == "propertypes":
                    status = "impropertypes"
                    bonfile2.write(f"{line}\n")
                    bonfile1.write(f"{line}\n")
                    bonfile1.write("[ dihedraltypes ] ; improper stuff\n")
                    bonfile1.write(";i   j   k   l	   func\n")
                elif status == "impropertypes":
                    status = "gap_after_improper"
                    bonfile1.write(f"{line}\n")
                elif status == "gap_after_improper":
                    pass
                elif status == "atomtypes":
                    status = "end"
                    nonfile.close()
                    docfile.write(f"{line}\n")
            elif "MOD4" in line and "RE" in line:
                status = "atomtypes"
                bonfile1.close()
                nonfile.write("[ atomtypes ]\n")
                nonfile.write("; name      at.num  mass     charge ptype  sigma      epsilon\n")
            else:
                if status == "initial":
                    status = "masstypes"
                    docfile.write(f"{line}\n")
                    atpfile.write(f"{line}\n")
                elif status == "masstypes":
                    parts = line.split()
                    mass[parts[0]] = float(parts[1])
                    len_line = len(line)
                    pos = line.find(parts[3])
                    if pos > 0:
                        line = f"{parts[0]:<3}{float(parts[1]):>10.5f} ; {line[pos:len_line]}"
                    else:
                        line = f"{parts[0]:<3}{float(parts[1]):>10.5f}\n"
                    atpfile.write(f"{line}\n")
                elif status == "hydrophilic":
                    status = "bondtypes"
                elif status == "bondtypes":
                    parts = line.replace("-", " ").split()
                    func = 1
                    pos = line.find(parts[4])
                    len_line = len(line)
                    if pos > 0:
                        line = f"{parts[0]:<10}{parts[1]:<10}{func:2}{0.1*float(parts[3]):>10.5f}{100*2*cal*float(parts[2]):>12.1f} ; {line[pos:len_line]}"
                    else:
                        line = f"{parts[0]:<10}{parts[1]:<10}{func:2}{0.1*float(parts[3]):>10.5f}{100*2*cal*float(parts[2]):>12.1f}\n"
                    bonfile1.write(f"{line}\n")
                elif status == "angletypes":
                    parts = line.replace("-", " ").split()
                    func = 1
                    pos = line.find(parts[5])
                    len_line = len(line)
                    if pos > 0:
                        line = f"{parts[0]:<10}{parts[1]:<10}{parts[2]:<10}{func:2}{float(parts[4]):>10.3f}{2*cal*float(parts[3]):>10.3f} ; {line[pos:len_line]}"
                    else:
                        line = f"{parts[0]:<10}{parts[1]:<10}{parts[2]:<10}{func:2}{float(parts[4]):>10.3f}{2*cal*float(parts[3]):>10.3f}\n"
                    bonfile1.write(f"{line}\n")
                elif status == "propertypes":
                    parts = line.replace("-", " ").split()
                    func = 9
                    try:
                        pos = line.find(parts[8])
                        len_line = len(line)
                        if float(parts[4]) <= 0:
                            print("error: IDIVF <= 0 at the line")
                            print(line)
                        if pos > 0:
                            line = f"{parts[0]:<10}{parts[1]:<10}{parts[2]:<10}{parts[3]:<10}{func:2}{float(parts[6]):>8.1f}{cal*float(parts[5])/float(parts[4]):>12.5f}{abs(float(parts[7])):>5.0f} ; {line[pos:len_line]}"
                        else:
                            line = f"{parts[0]:<10}{parts[1]:<10}{parts[2]:<10}{parts[3]:<10}{func:2}{float(parts[6]):>8.1f}{cal*float(parts[5])/float(parts[4]):>12.5f}{abs(float(parts[7])):>5.0f}\n"
                    except IndexError: 
                        line = f"{parts[0]:<10}{parts[1]:<10}{parts[2]:<10}{parts[3]:<10}{func:2}{float(parts[6]):>8.1f}{cal*float(parts[5])/float(parts[4]):>12.5f}{abs(float(parts[7])):>5.0f}\n"
                    bonfile2.write(f"{line}\n")
                elif status == "impropertypes":
                    parts = line.replace("-", " ").split()
                    func = 4
                    try:
                        pos = line.find(parts[7])
                        len_line = len(line)
                        if pos > 0:
                            line = f"{parts[0]:<10}{parts[1]:<10}{parts[2]:<10}{parts[3]:<10}{func:2}{float(parts[5]):>8.1f}{cal*float(parts[4]):>12.5f}{float(parts[6]):>5.0f} ; {line[pos:len_line]}"
                        else:
                            line = f"{parts[0]:<10}{parts[1]:<10}{parts[2]:<10}{parts[3]:<10}{func:2}{float(parts[5]):>8.1f}{cal*float(parts[4]):>12.5f}{float(parts[6]):>5.0f}\n"
                    except IndexError:
                        line = f"{parts[0]:<10}{parts[1]:<10}{parts[2]:<10}{parts[3]:<10}{func:2}{float(parts[5]):>8.1f}{cal*float(parts[4]):>12.5f}{float(parts[6]):>5.0f}\n"
                    bonfile1.write(f"{line}\n")
                elif status == "gap_after_improper":
                    pass
                elif status == "atomtypes":
                    parts = line.split()
                    atnum = 0
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
                    line = f"{parts[0]:<10}{atnum:>3}{mass[parts[0]]:>8.3f}  0.0  A {float(parts[1])*0.1*(2**(5/6)):.5e} {cal*float(parts[2]):.5e} ; {line[41:]}"
                    nonfile.write(f"{line}\n")
                elif status == "end":
                    if parts[0][:3] != "END":
                        len_line = len(line)
                        line = f"; {line[:len_line]}"
                        docfile.write(f"{line}\n")
