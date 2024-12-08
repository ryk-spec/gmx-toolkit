import os
from gmx_toolkit.runner.gmx_run import run_command_with_input
import opt_axes
import pandas as pd
import mdtraj as md
import numpy as np

def xvg_to_df(xvg):

    with open(xvg) as f:
        lines = f.readlines()
    lines = [line.split() for line in lines]
    df = pd.DataFrame(lines,columns=["X","Y"])
    return df

class GmxEnergy:

    def __init__(self,edr):

        self.edr = edr 
        self.xvg = None

    
    def run(self,thermo="density"):

        self.xvg = os.path.splitext(self.edr)[0]+f"_{thermo}.xvg"
        command = f"gmx energy -f {self.edr} -o {self.xvg} -xvg none"
        run_command_with_input(
            command=command,
            input_data="density"
        )
        self.df = xvg_to_df(self.xvg)
    
    
class WaxsCalculator:

    def __init__(self,gro,trr):

        self.gro = gro
        self.trr = trr

    def assign_asf(self):

        asf_a = {}
        asf_b = {}
        asf_c = {}
        asf_a["H"] = np.array([0.493002,0.322912,0.140191,0.04081])
        asf_b["H"] = np.array([10.5109,26.1257,3.14236,57.7997])
        asf_c["H"] = np.array([0.003038])
        asf_a["C"] = np.array([2.31,1.02,1.5886,0.865])
        asf_b["C"] = np.array([20.8439,10.2075,0.5687,51.6512])
        asf_c["C"] = np.array([0.2156])
        data = [12.2126, 0.0057, 3.1322, 9.8933, 2.0125, 28.9975, 1.1663, 0.5826, -11.529]
        # Reshape into arrays for a, b, c
        asf_a["N"] = np.array(data[:-1:2])  # Elements at even indices
        asf_b["N"] = np.array(data[1::2])  # Elements at odd indices except the last
        asf_c["N"] = np.array([data[-1]])  # Last element
        data = [3.0485, 13.2771, 2.2868, 5.7011, 1.5463, 0.3239, 0.867, 32.9089, 0.2508]
        asf_a["O"] = np.array(data[:-1:2])  # Elements at even indices
        asf_b["O"] = np.array(data[1::2])  # Elements at odd indices except the last
        asf_c["O"] = np.array([data[-1]])  # Last element
        trajectory = md.load(self.gro)  # 必要なファイル名に変更
        topology = trajectory.topology
        # 条件に基づく原子番号を取得
        # 例: 元素が 'C' の原子番号を取得
        self.asf_a_arr = np.stack([asf_a[atom.element.symbol] for atom in topology.atoms if atom.element.symbol != "VS"])
        self.asf_b_arr = np.stack([asf_b[atom.element.symbol] for atom in topology.atoms if atom.element.symbol != "VS"])
        self.asf_c_arr = np.stack([asf_c[atom.element.symbol] for atom in topology.atoms if atom.element.symbol != "VS"])
        print(
            self.asf_a_arr.shape,
            self.asf_b_arr.shape,
            self.asf_c_arr.shape
        )
        
