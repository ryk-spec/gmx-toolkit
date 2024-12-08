import os
from gmx_toolkit.runner.gmx_run import run_command_with_input
import opt_axes
import pandas as pd
import mdtraj as md
import numpy as np
import torch
from torch.utils.data import DataLoader, TensorDataset
import time
from scipy.interpolate import interp1d
import pandas as pd

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

    def assign_asf(self,atoms_idx):

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
        self.asf_a_arr = np.stack([asf_a[atom.element.symbol] for idx,atom in enumerate(topology.atoms) if idx in atoms_idx])
        self.asf_b_arr = np.stack([asf_b[atom.element.symbol] for idx,atom in enumerate(topology.atoms) if idx in atoms_idx])
        self.asf_c_arr = np.stack([asf_c[atom.element.symbol] for idx,atom in enumerate(topology.atoms) if idx in atoms_idx])

    def get_water_idx(self,topology : md.Topology):

        return topology.select("(symbol != VS) and water")

    def get_MOL_idx(self,topology : md.Topology):

        return topology.select("(symbol != VS) and (resname == MOL)")
        
    def get_all_idx(self,topology : md.Topology):

        return topology.select("(symbol != VS) and all")


    def mainroop(self,name,query="water",maxq=35,stride=10,end=100):
        start = time.time()
        frame : md.Trajectory
        print("begin to calculation")
        count = 0
        Iq_save_list = []
        for frame  in md.iterload(self.trr,top=self.gro,chunk=1,stride=stride):
            count += 1
            print(f"Frame number {frame.time[0]}",flush=True)
            if frame.time[0] >= end:
                break
            topology = frame.topology
            atoms_idx = None
            if query == "water":
                atoms_idx = self.get_water_idx(topology=topology)
            elif query == "MOL":
                atoms_idx = self.get_MOL_idx(topology=topology)
            elif query == "all":
                atoms_idx = self.get_all_idx(topology=topology)
            coords = frame.xyz[0,atoms_idx,:]
            L = frame.unitcell_lengths[0,0]
            maxq_int = int(maxq/(2*np.pi/L))
            qs = np.array([[i,j,k] 
                  for i in np.arange(-maxq_int,maxq_int+1) 
                  for j in np.arange(-maxq_int,maxq_int+1) 
                  for k in np.arange(-maxq_int,maxq_int+1)],dtype=np.float32)
            qs *= 2*np.pi/L
            self.assign_asf(atoms_idx)
            F=self.calc_F(qs,coords)
            q_norm = np.linalg.norm(qs,axis=1)
            q_id = q_norm / (2*np.pi/L)
            Iqs = []
            for i in range(1,maxq_int):
                F_selected = F[(q_id>=i-0.5) & (q_id<i+0.5)]
                Iq = np.mean(np.abs(F_selected)**2)
                Iq /= (L*10**(-7))**3
                Iq *= (2.818e-13)**2
                Iqs.append(Iq)
                print(i*2*np.pi/L,Iq)
            Iq = np.array(Iqs)
            waxs_interp = interp1d(
                x=2*np.pi/L*np.arange(1,maxq_int),
                y=Iq,kind="cubic",fill_value="extrapolate"
            )
            Iq_save_list.append(
                waxs_interp(np.linspace(0,maxq,100))
            )
        Iq = np.stack(Iq_save_list,axis=-1)
        Iq = np.mean(Iq,axis=1)
        print(np.stack([np.linspace(0,maxq,100),Iq],axis=1))
        df = pd.DataFrame(
            np.stack([np.linspace(0,maxq,100),Iq],axis=1),
            columns=["x","y"]
        )
        df = df.map(lambda x: '{0:.5f}'.format(x))
        df.to_csv(
            os.path.splitext(self.trr)[0]+"_"+name+".xvg",
            sep=" ",columns=None,index=None
        )
        end = time.time()
        print(f"\nend {end-start:.2f} second")

        
    def calc_F(self,qs,coords):
        a = torch.from_numpy(self.asf_a_arr).to(torch.float32)  # a_{j,k}
        b = torch.from_numpy(self.asf_b_arr).to(torch.float32)  # (N,4)
        c = torch.from_numpy(self.asf_c_arr).to(torch.float32)  # a_{j,k}
        r = torch.from_numpy(coords).to(torch.float32)  # r_j (位置ベクトル)
        q = torch.from_numpy(qs).to(torch.float32)  # r_j (位置ベクトル)
        dataset = TensorDataset(q)
        dataset2 = TensorDataset(a,b,c,r)
        # DataLoader を使用して、データをバッチごとにロード
        batch_size =2**13
        dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=False)
        dataloader2 = DataLoader(dataset2, batch_size=batch_size, shuffle=False)
        F_real_list = []  # F_real を保存するリスト
        F_imag_list = []
        r_gpu = r.to("cuda")
        for batch_q in dataloader:
            batch_q = batch_q[0]
            batch_q = batch_q.to(device="cuda")
            q_num = batch_q.size()[0]
            # e^{-||q||} の計算
            q_norm = torch.norm(batch_q, dim=1) # (q_num)
            q_norm = q_norm.unsqueeze(0)
            inner_sum_list = []
            for batch_a,batch_b,batch_c,batch_r in dataloader2:
                batch_a = batch_a.to(device="cuda")
                batch_b = batch_b.to(device="cuda")
                batch_c = batch_c.to(device="cuda")
                batch_r = batch_r.to(device="cuda")
                b_extend = batch_b.unsqueeze(-1)
                exp_q_norm = torch.exp(-torch.matmul(b_extend,q_norm.pow(2) / (16 * torch.pi**2 * 100)))  # (N,4,q_num)
                # jごとの内部和を計算
                inner_sum = torch.sum(torch.stack([batch_a for _ in range(q_num)],dim=-1) * exp_q_norm, dim=1) + torch.cat([batch_c for _ in range(q_num)],dim=-1)  # (N,q_num)
                inner_sum_list.append(inner_sum)
            inner_sum = torch.cat(inner_sum_list, dim=0)
            # q ⋅ r_j の計算
            phase = torch.matmul(r_gpu, batch_q.t())  # qとrの内積 (num_atoms, num_q)
            phase_real = torch.cos(phase)
            phase_imag = torch.sin(phase)
            # F の実部と虚部を更新
            F_real = torch.sum(inner_sum * phase_real, dim=0)  # (num_atoms,)と(num_atoms, num_q)の積
            F_imag = torch.sum(inner_sum * phase_imag, dim=0)  # (num_atoms,)と(num_atoms, num_q)の積
            F_real_list.append(F_real)
            F_imag_list.append(F_imag)

        del q_norm, exp_q_norm, b_extend,inner_sum
        del phase, phase_real, phase_imag
        torch.cuda.empty_cache()
        # 結果の複素数形式
        F_real = torch.cat(F_real_list, dim=0)  # (num_q,)
        F_imag = torch.cat(F_imag_list, dim=0)  # (num_q,)
        F = F_real + 1j * F_imag
        F = F.to("cpu").numpy()
        del a, b, c, r
        torch.cuda.empty_cache()
        # 結果を表示
        return F
