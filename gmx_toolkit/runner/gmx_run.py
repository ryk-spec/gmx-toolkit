import os
import subprocess
import numpy as np
import mdtraj as md
import re
import scipy.constants as const

def run_command_with_input(command, input_data=None):
    try:
        command = command.split()
        # コマンドを実行し、標準入力を与える
        result = subprocess.run(
            command,  # 実行したいコマンドと引数
            check=True,  # コマンドが失敗した場合に例外を発生させる
            stdout=subprocess.PIPE,  # 標準出力をキャプチャ
            stderr=subprocess.PIPE,  # 標準エラー出力をキャプチャ
            text=True,  # 出力を文字列として取得
            input=input_data  # 標準入力データを渡す
        )
    except subprocess.CalledProcessError as e:
        # エラーが発生した場合、エラーメッセージを表示
        print(f"{e.stderr}")  # 標準エラー出力
        raise  # 例外を再度発生させてPython側でもエラーを起こす
    else:
        # 成功した場合の処理（標準出力を使いたい場合）
        return result.stdout  # 成功した場合に標準出力を返す

class ForceFieldSetup:
    """
    分子シミュレーション用の力場セットアップとトポロジー作成を管理するクラス。
    
    属性:
        ffdir (str): 力場ファイルが格納されているディレクトリ。
        data_dict (dict): パースした .itp ファイルのデータをセクションごとに整理した辞書。
        mole_name (str): .itp ファイルから抽出した分子名。
        itp_path (str): 使用中の .itp ファイルのパス。
    """

    def __init__(self, ffdir=None):
        """
        ForceFieldSetup クラスを初期化する。
        
        引数:
            ffdir (str, オプション): 力場ディレクトリのパス。指定がない場合は、現在のディレクトリ内の "gaff.ff" を使用。
        """
        # ディレクトリが指定されていない場合は、デフォルトでリポジトリ内の "gaff.ff" を使用
        if ffdir is None:
            ffdir = os.path.join(os.path.dirname(__file__), "gaff.ff")
        self.ffdir = ffdir
        self.data_dict = None
        self.mole_name = None
        self.itp_path = None

    def itp_parser(self, itp_path):
        """
        .itp ファイルをパースして、その内容を辞書に整理する。
        
        引数:
            itp_path (str): パースする .itp ファイルのパス。
            
        例外:
            FileNotFoundError: 指定した .itp ファイルが存在しない場合に発生。
        """
        self.itp_path = itp_path
        data_dict = {}

        # .itp ファイルを読み込む
        with open(itp_path, 'r', encoding='utf-8') as file:
            content = file.read()

        # セクションごとにコンテンツを分割（例: [ moleculetype ], [ atoms ]）
        sections = re.split(r'(\[.*?\])', content)
        # 最初の要素が空の場合は削除
        if "[" not in sections[0].strip() :
            sections.pop(0)

        # セクション名とその内容をペアにして辞書に格納
        for i in range(0, len(sections), 2):
            key = sections[i].strip()
            value = sections[i + 1].strip() if i + 1 < len(sections) else ''

            # 重複するキーがあれば、キーに番号を付けてユニークにする
            original_key = key
            counter = 2
            while key in data_dict:
                key = f"{original_key}_{counter}"
                counter += 1

            data_dict[key] = value

        self.data_dict = data_dict

    def get_moleculename(self):
        """
        パースした .itp ファイルから分子名を抽出して保存する。
        
        例外:
            AssertionError: `self.data_dict` が None の場合、
            つまり .itp ファイルがまだパースされていない場合に発生。
        """
        assert self.data_dict is not None, "\
            最初に `itp_parser` を使って .itp ファイルをパースしてください。"
        # [ moleculetype ] セクションから分子名を抽出
        # print(self.data_dict.keys())
        self.mole_name = self.data_dict["[ moleculetype ]"]
        self.mole_name = self.mole_name.split("\n")[-1].split()[0]

    def write_itp(self, data_dict, keys, itp_path):
        """
        修正した .itp ファイルを保存する。指定したキーに一致するセクションは除外される。
        
        引数:
            data_dict (dict): 修正された .itp のデータ。
            keys (list of str): 除外するセクションのキーワードのリスト。
            itp_path (str): 修正した .itp ファイルを保存するパス。
        """
        with open(itp_path, mode="w", encoding="utf-8") as f:
            to_remove = set()
            # 除外すべきセクションを特定
            for query in keys:
                for key in data_dict.keys():
                    if query in key:
                        to_remove.add(key)
            # 残りのセクションをファイルに書き込む
            for key, value in data_dict.items():
                lines = value.split("\n")
                if key not in to_remove:
                    if key == "[ bonds ]":
                        for idx,line in enumerate(lines):
                            if len(line) == 0:
                                continue
                            if line[0] == ";":
                                continue
                            else:
                                lines[idx] = " "+"\t".join(line.split()[:2])
                    if key == "[ angles ]":
                        for idx,line in enumerate(lines):
                            if len(line) == 0:
                                continue
                            if line[0] == ";":
                                continue
                            else:
                                lines[idx] = " "+"\t".join(line.split()[:3])
                    if key == "[ dihedrals ]":
                        for idx,line in enumerate(lines):
                            if len(line) == 0:
                                continue
                            if line[0] == ";":
                                continue
                            else:
                                lines[idx] = " "+"\t".join(line.split()[:5])
                    lines.append(" ")
                    value = "\n".join(lines)
                    f.write(key + "\n")
                    f.write(value + " \n")

    @classmethod
    def write_top(cls, top_path, itp_paths, nums, molecules, name_of_run):
        """
        トポロジー (.top) ファイルを生成し保存する。
        
        引数:
            top_path (str): 保存する .top ファイルのパス。
            itp_paths (list of str): インクルードする .itp ファイルのパスのリスト。
            nums (list of int): 各分子の分子数。
            molecules (list of str): 分子の名前。
            name_of_run (str): シミュレーションシステムの名前。
        """
        contents = []
        # デフォルトセクションを追加
        # contents.append("[ defaults ]")
        # contents.append("1 2 yes 0.5 0.8333")
        
        # .itp ファイルをインクルード
        for itp in itp_paths:
            contents.append(f'#include "{itp}"')
        
        # システム名を追加
        contents.append('[ system ]')
        contents.append(f'{name_of_run}')
        
        # 分子セクションを追加
        contents.append('[ molecules ]')
        for molecule, num in zip(molecules, nums):
            if num is not None:
                contents.append(f'{molecule} {num}')
        
        # .top ファイルを書き込む
        with open(top_path, "w", encoding="utf-8") as f:
            f.write("\n".join(contents) + "\n")



class InitialStructure:
    """
    分子シミュレーションおよび溶媒和設定を管理するクラス。

    属性:
    ----------
    density : float
        システムの密度 (g/cm³)。

    メソッド:
    -------
    get_M(gro_path, resid)
        指定された .gro ファイルの最初の分子の分子量を計算します。
    register_solv(solv_pdb)
        PDBファイル内の残基名を調整し、溶媒和に適した形式に変換します。
    solv_setup(solv_gro_path, solv_itp_path, box_size=5)
        溶媒和パラメータを設定し、GROMACSを使用して溶媒和システムを生成します。
    """

    def __init__(self, density):
        """
        InitialStructureオブジェクトを密度を指定して初期化します。

        パラメータ:
        ----------
        density : float
            システムの密度 (g/cm³)。
        """
        self.density = density

    @staticmethod
    def get_M(gro_path, resid):
        """
        指定された .gro ファイル内の特定の分子の分子量を計算します。

        パラメータ:
        ----------
        gro_path : str
            分子情報を含む .gro ファイルへのパス。
        resid : int
            分子量を計算する対象の残基インデックス。

        戻り値:
        -------
        float
            分子量（原子質量単位：amu）。
        """
        # .groファイルからトラジェクトリを読み込む
        trajectory: md.Trajectory = md.load(gro_path)
        top: md.Topology = trajectory.topology
        
        # 指定された残基の原子質量を取得
        first_molecule_masses = [atom.element.mass for atom in top.atoms 
            if atom.residue.index == resid]
        
        # 原子質量の合計を分子量として計算
        molecular_weight = sum(first_molecule_masses)
        return molecular_weight

    def register_solv(self, solv_pdb):
        """
        PDBファイル内の残基名を修正し、溶媒和の準備を行います。

        パラメータ:
        ----------
        solv_pdb : str
            溶媒のPDBファイルへのパス。
        """
        # PDBファイルを読み込む
        traj = md.load(solv_pdb)
        
        # 全ての残基名を "SOL" に変更
        for residue in traj.topology.residues:
            residue.name = "SOL"
        
        # 修正後のPDBファイルを保存
        output_pdb = os.path.splitext(solv_pdb)[0] + "_registered.pdb"
        traj.save(output_pdb)
        
        # GROMACSのinsert-moleculesコマンドを生成
        insert = f"gmx insert-molecules -box 5 5 5 \
            --ci {output_pdb} -nmol 4000 -try 3 \
            -o {os.path.splitext(solv_pdb)[0]}.gro"
        
        # コマンドを実行
        run_command_with_input(insert)

    def solv_setup(self, solv_gro_path, solv_itp_path, box_size=5):
        """
        GROMACSツールを使用して溶媒和を設定します。

        パラメータ:
        ----------
        solv_gro_path : str
            溶媒の .gro ファイルへのパス。
        solv_itp_path : str
            溶媒の .itp ファイルへのパス。
        box_size : float, optional
            シミュレーションボックスのサイズ (デフォルトは5 nm)。
        """
        # .groファイルの最初の分子の分子量を計算
        M_0 = InitialStructure.get_M(solv_gro_path, 0)
        
        # ボックスサイズをnm³からcm³に変換し、必要な分子数を計算
        N_0_float = (box_size * 1e-7)**3  # 体積 (cm³)
        N_0_float *= self.density  # 質量 (g)
        N_0_float *= const.Avogadro  # モルに変換
        N_0_float /= M_0  # 分子数を計算
        N_0 = int(N_0_float)  # 整数に丸める
        
        # 力場をセットアップ
        solv = ForceFieldSetup(ffdir=None)
        solv.itp_parser(itp_path=solv_itp_path)
        solv.get_moleculename()
        
        # 溶媒和システムのトポロジーファイルを作成
        ForceFieldSetup.write_top(
            top_path=f"{os.path.basename(os.path.splitext(solv_itp_path)[0])}_solvated.top",
            itp_paths=[solv.ffdir + "/forcefield.itp", solv.itp_path],
            nums=[None],
            molecules=[solv.mole_name],
            name_of_run=f"{os.path.basename(os.path.splitext(solv_itp_path)[0])}"
        )
        
        # GROMACSのsolvateコマンドを生成
        solvate = f"gmx solvate -box {box_size} {box_size} {box_size} \
            -cs {solv_gro_path} \
            -o {os.path.basename(os.path.splitext(solv_itp_path)[0])}_solvated.gro \
            -p {os.path.basename(os.path.splitext(solv_itp_path)[0])}_solvated.top -maxsol {N_0}"
        
        # コマンドを実行
        run_command_with_input(solvate)

    def soln_setup(self, solute_gro_path,solute_itp_path,solute_itp_path_out,solv_gro_path, solv_itp_path, wf,box_size=5):
        """
        GROMACSツールを使用して溶媒和を設定します。

        パラメータ:
        ----------
        solute_gro_path : str
            溶質の .gro ファイルへのパス。
        solute_itp_path : str
            溶質の .itp ファイルへのパス。
        solute_itp_path : str
            溶質の .itp ファイルへのパス。(atomtypesを省き出力)
        solv_gro_path : str
            溶媒の .gro ファイルへのパス。
        solv_itp_path : str
            溶媒の .itp ファイルへのパス。
        box_size : float, optional
            シミュレーションボックスのサイズ (デフォルトは5 nm)。
        """
        # .groファイルの最初の分子の分子量を計算
        M_0 = InitialStructure.get_M(solv_gro_path, 0)
        M_1 = InitialStructure.get_M(solute_gro_path, 0)
        # ボックスサイズをnm³からcm³に変換し、必要な分子数を計算
        N_1_float = (box_size * 1e-7)**3  # 体積 (cm³)
        N_1_float *= self.density  # 質量 (g)
        N_1_float *= wf
        N_1_float *= const.Avogadro  # モルに変換
        N_1_float /= M_1  # 分子数を計算
        N_1 = int(N_1_float)  # 整数に丸める
        N_0 = N_1 * M_1 / wf * (1-wf) / M_0
        N_0 = int(N_0)
        # 力場をセットアップ
        solv = ForceFieldSetup(ffdir=None)
        solv.itp_parser(itp_path=solv_itp_path)
        solv.get_moleculename()
        solute = ForceFieldSetup(ffdir=None)
        solute.itp_parser(itp_path=solute_itp_path)
        solute.write_itp(
            data_dict=solute.data_dict,
            keys=["atomtypes"],
            itp_path=solute_itp_path_out
        )
        solute.get_moleculename()
        # 溶媒和システムのトポロジーファイルを作成
        ForceFieldSetup.write_top(
            top_path=f"{os.path.basename(os.path.splitext(solute_itp_path_out)[0])}_{os.path.basename(os.path.splitext(solv_itp_path)[0])}_solvated.top",
            itp_paths=[solv.ffdir + "/forcefield.itp",solute_itp_path_out,solv.itp_path],
            nums=[N_1,None],
            molecules=[solute.mole_name,solv.mole_name],
            name_of_run=f"{os.path.basename(os.path.splitext(solute_itp_path_out)[0])} in {os.path.basename(os.path.splitext(solv_itp_path)[0])}"
        )
        
        # GROMACSのsolvateコマンドを生成
        insert = f"gmx insert-molecules -box {box_size} {box_size} {box_size} \
            -ci {solute_gro_path} -nmol {N_1} \
            -o {os.path.basename(os.path.splitext(solute_itp_path_out)[0])}_randomized.gro"
        solvate = f"gmx solvate -cp {os.path.basename(os.path.splitext(solute_itp_path_out)[0])}_randomized.gro \
            -cs {solv_gro_path} \
            -o {os.path.basename(os.path.splitext(solute_itp_path_out)[0])}_{os.path.basename(os.path.splitext(solv_itp_path)[0])}_solvated.gro \
            -p {os.path.basename(os.path.splitext(solute_itp_path_out)[0])}_{os.path.basename(os.path.splitext(solv_itp_path)[0])}_solvated.top\
            -maxsol {N_0}"
        
        # コマンドを実行
        run_command_with_input(insert)
        run_command_with_input(solvate)


class RunThrow:
    """
    GROMACSシミュレーションの準備と実行を行うクラス。

    このクラスは、エネルギー最小化（EM）、NVT、NPT、製造シミュレーションの各ステージを実行するために使用されます。
    また、各シミュレーションに必要な.mdpファイルを生成するためのメソッドも提供しています。
    """

    def __init__(self, topfile, init_grofile):
        """
        RunThrowクラスのインスタンスを初期化します。

        :param topfile: トポロジーファイル（.top）のパス（str）
        :param init_grofile: 初期座標ファイル（.gro）のパス（str）
        """
        self.topfile = topfile
        self.init_grofile = init_grofile

    def run_EM(self):
        """
        エネルギー最小化（EM）シミュレーションを実行します。

        GROMACSのエネルギー最小化（EM）を行い、生成される最小化された構造を使用して次のシミュレーションを準備します。
        """
        RunThrow.generate_em_mdp()
        tpr = os.path.splitext(self.topfile)[0]
        tpr = tpr.replace("solvated", "EM") + ".tpr"
        grompp = f"gmx grompp -c {self.init_grofile} -p {self.topfile} -f em.mdp -o {tpr}"
        mdrun = f"gmx mdrun -s {tpr} -deffnm {os.path.basename(tpr)}"
        run_command_with_input(grompp)
        run_command_with_input(mdrun)

    def run_nvt(self, T, run_ns=0.5,grps=["SOL"]):
        """
        NVT（定温定容）シミュレーションを実行します。

        GROMACSのNVTシミュレーションを実行し、指定された温度でシステムを熱化します。

        :param T: 目標温度（K）
        :param run_ns: シミュレーション時間（ns）
        """
        RunThrow.generate_nvt_mdp(temperature=T, run_ns=run_ns,grps=grps)
        tpr = os.path.splitext(self.topfile)[0]
        T = str(T).replace(".", "p")
        prev_tpr = tpr.replace("solvated", "EM")
        tpr = tpr.replace("solvated", f"{T}_nvt")
        grompp = f"gmx grompp -c {prev_tpr}.tpr.gro -p {self.topfile} -f nvt.mdp -o {tpr}.tpr"
        mdrun = f"gmx mdrun -s {tpr}.tpr -deffnm {os.path.basename(tpr)}.tpr -cpi {tpr}.tpr.cpt"
        run_command_with_input(grompp)
        run_command_with_input(mdrun)

    def run_npt(self, T, run_ns=5,grps=["SOL"]):
        """
        NPT（定温定圧）シミュレーションを実行します。

        GROMACSのNPTシミュレーションを実行し、指定された温度と圧力でシステムを緩和します。

        :param T: 目標温度（K）
        :param run_ns: シミュレーション時間（ns）
        """
        RunThrow.generate_npt_mdp(temperature=T, run_ns=run_ns,grps=grps)
        tpr = os.path.splitext(self.topfile)[0]
        T = str(T).replace(".", "p")
        prev_tpr = tpr.replace("solvated", f"{T}_nvt")
        tpr = tpr.replace("solvated", f"{T}_npt")
        grompp = f"gmx grompp -c {prev_tpr}.tpr.gro -p {self.topfile} -f npt.mdp -o {tpr}.tpr"
        mdrun = f"gmx mdrun -s {tpr}.tpr -deffnm {os.path.basename(tpr)}.tpr -cpi {tpr}.tpr.cpt"
        run_command_with_input(grompp)
        run_command_with_input(mdrun)

    def run_production(self, T, run_ns=5,grps=["SOL"]):
        """
        生産シミュレーションを実行します。

        GROMACSの生産シミュレーションを実行し、指定された温度でシステムを長時間シミュレートします。

        :param T: 目標温度（K）
        :param run_ns: シミュレーション時間（ns）
        """
        RunThrow.generate_production_mdp(temperature=T, run_ns=run_ns,grps=grps)
        tpr = os.path.splitext(self.topfile)[0]
        T = str(T).replace(".", "p")
        prev_tpr = tpr.replace("solvated", f"{T}_npt")
        tpr = tpr.replace("solvated", f"{T}_production")
        grompp = f"gmx grompp -c {prev_tpr}.tpr.gro -p {self.topfile} -f production.mdp -o {tpr}.tpr"
        mdrun = f"gmx mdrun -s {tpr}.tpr -deffnm {os.path.basename(tpr)}.tpr -cpi {tpr}.tpr.cpt"
        run_command_with_input(grompp)
        run_command_with_input(mdrun)

    @staticmethod
    def write_mdp(filename, parameters):
        """
        汎用的な.mdpファイルを書き込む関数。

        :param filename: 出力する.mdpファイルの名前（str）
        :param parameters: .mdpに書き込むパラメータ（dict）
        """
        with open(filename, 'w', encoding='utf-8') as f:
            for key, value in parameters.items():
                f.write(f"{key} = {value}\n")
        print(f"{filename} を作成しました。")

    @staticmethod
    def generate_em_mdp(output_file="em.mdp"):
        """
        エネルギー最小化用の.mdpファイルを作成します。

        :param output_file: 出力ファイル名（str）
        """
        em_params = {
            "integrator": "steep",        # 最急降下法
            "nsteps": 50000,             # 最大ステップ数
            "emtol": 10.0,             # エネルギー収束基準
            "emstep": 0.01,              # 最小化ステップサイズ
            "cutoff-scheme": "Verlet",
            "nstlist": 10,
            "coulombtype": "PME",
            "rcoulomb": 1.5,
            "rvdw": 1.5
        }
        RunThrow.write_mdp(output_file, em_params)

    @staticmethod
    def generate_nvt_mdp(output_file="nvt.mdp", temperature=300, run_ns=0.5, grps=['SOL']):
        """
        NVT緩和用の.mdpファイルを作成します。

        :param output_file: 出力ファイル名（str）
        :param temperature: 目標温度（K）
        :param run_ns: シミュレーション時間（ns）
        :param grps: 速度生成の対象となるグループ（list）
        """
        nsteps = run_ns * 1000 / 0.0005
        nsteps = int(nsteps)
        nvt_params = {
            "integrator": "md",          # 分子動力学
            "nsteps": nsteps,            # シミュレーションステップ数
            "dt": 0.0005,                 # タイムステップ（ps）
            "cutoff-scheme": "Verlet",
            "nstlist": 10,
            "coulombtype": "PME",
            "rcoulomb": 1.0,
            "rvdw": 1.0,
            "tcoupl": "V-rescale",       # 温度カップリング法
            "pcoupl": "no",              # 圧力カップリングなし（NVT）
            "gen_vel": "yes",            # 初期速度を生成
            "gen_temp": temperature,
            "gen_seed": -1,
        }
        nvt_params["tc-grps"] = "\t".join(grps)
        nvt_params["tau-t"] = "\t".join(["0.2" for _ in grps])
        nvt_params["ref-t"] = "\t".join([f"{temperature}" for _ in grps])
        RunThrow.write_mdp(output_file, nvt_params)

    @staticmethod
    def generate_npt_mdp(output_file="npt.mdp", temperature=300, run_ns=5, pressure=1.0, grps=["SOL"]):
        """
        NPT緩和用の.mdpファイルを作成します。

        :param output_file: 出力ファイル名（str）
        :param temperature: 目標温度（K）
        :param pressure: 目標圧力（bar）
        :param run_ns: シミュレーション時間（ns）
        :param grps: 速度生成の対象となるグループ（list）
        """
        nsteps = run_ns * 1000 / 0.002
        nsteps = int(nsteps)
        npt_params = {
            "integrator": "md",
            "nsteps": nsteps,
            "dt": 0.002,
            "cutoff-scheme": "Verlet",
            "nstlist": 10,
            "coulombtype": "PME",
            "rcoulomb": 1.5,
            "rvdw": 1.5,
            "tcoupl": "V-rescale",
            "pcoupl": "C-rescale",  # 圧力カップリング
            "pcoupltype": "isotropic",
            "tau_p": 2.0,
            "ref_p": pressure,
            "compressibility": 4.5e-5,
            "gen_vel": "no",
            "constraints": "all-bonds"
        }
        npt_params["tc-grps"] = "\t".join(grps)
        npt_params["tau-t"] = "\t".join(["0.5" for _ in grps])
        npt_params["ref-t"] = "\t".join([f"{temperature}" for _ in grps])
        RunThrow.write_mdp(output_file, npt_params)

    @staticmethod
    def generate_production_mdp(output_file="production.mdp", run_ns=5, temperature=300, grps=['SOL']):
        """
        生産シミュレーション用の.mdpファイルを作成します。

        :param output_file: 出力ファイル名（str）
        :param run_ns: シミュレーション時間（ns）
        :param temperature: 目標温度（K）
        :param grps: 速度生成の対象となるグループ（list）
        """
        nsteps = run_ns * 1000 / 0.002
        nsteps = int(nsteps)
        production_params = {
            "integrator": "md",
            "nsteps": nsteps,
            "dt": 0.002,
            "cutoff-scheme": "Verlet",
            "nstlist": 10,
            "coulombtype": "PME",
            "rcoulomb": 1.0,
            "rvdw": 1.0,
            "tcoupl": "nose-hoover",
            "pcoupl": "Parrinello-Rahman",
            "pcoupltype": "isotropic",
            "tau_p": 5.0,
            "ref_p": 1.0,
            "compressibility": 4.5e-5,
            "gen_vel": "no",
            "nstxout": 1000,
            "nstvout": 1000,
            "nstenergy": 1000,
            "nstlog": 1000,
            "constraints": "all-bonds"
        }
        production_params["tc-grps"] = "\t".join(grps)
        production_params["tau-t"] = "\t".join(["0.5" for _ in grps])
        production_params["ref-t"] = "\t".join([f"{temperature}" for _ in grps])
        RunThrow.write_mdp(output_file, production_params)

def solv_MD(solv_gro,solv_itp,box_size=5,nvt_ns=0.5,npt_ns=5,production_ns=5,T=301.15,solv_pdb="None"):
    """
    溶媒のMDを回すメソッド

    Args:
        solv_gro (str): 溶媒和された .gro
        solv_itp (str): 溶媒の .itp
        solv_pdb (str, optional): 溶媒1分子の .pdb 指定するとsolv_groを生成してくれる. Default は None.
        box_size (float): 箱のサイズ nm
        nvt_ns (float) : NVT緩和 ns
        npt_ns (float) : NPT緩和 ns
        production_ns (float) : production run ns
        T (float): 絶対温度
    """
    solv = InitialStructure(0.9)
    if solv_pdb != "None":
        solv.register_solv(
            solv_pdb=solv_pdb
        )
        solv_gro = os.path.splitext(solv_pdb)[0]+".pdb"
    solv.solv_setup(
        solv_gro_path=solv_gro,
        solv_itp_path=solv_itp,
        box_size=box_size
    )
    topfile = os.path.basename(solv_itp)
    topfile = os.path.splitext(topfile)[0]+"_solvated.top"
    grofile = os.path.basename(solv_itp)
    grofile = os.path.splitext(grofile)[0]+"_solvated.gro"
    solv_run = RunThrow(
        topfile=topfile,
        init_grofile=grofile
    )
    solv_run.run_EM()
    solv_run.run_nvt(T=T,run_ns=nvt_ns)
    solv_run.run_npt(T=T,run_ns=npt_ns)
    solv_run.run_production(T=T,run_ns=production_ns)

def soln_MD(solute_itp,solute_itp_out,solute_gro,solv_gro,solv_itp,box_size=8,wf=0.1,
            nvt_ns=0.5,npt_ns=10,production_ns=5,T=301.15,solv_pdb="None"):
    """
    溶液のMDを回すメソッド

    Args:
        solute_itp (_type_): 溶質の .itp
        solute_itp_out (_type_): 溶媒の .itp [ atomtypes ]を除去するのに用いる．ffnonbonded.itpにデータがないと怒られるよ！
        solute_gro (_type_): _description_
        solv_gro (str): 溶媒和された .gro
        solv_itp (str): 溶媒の .itp
        box_size (float): 箱のサイズ nm
        T (float): 絶対温度
        wf (float, optional): 溶質重量分率. Defaults to 0.1.
        nvt_ns (float) : NVT緩和 ns
        npt_ns (float) : NPT緩和 ns
        production_ns (float) : production run ns
        T (float, optional): _description_. Defaults to 301.15.
        solv_pdb (str, optional): 溶媒1分子の .pdb 指定するとsolv_groを生成してくれる. Default は None.
    """

    soln = InitialStructure(0.9)
    if solv_pdb != "None":
        soln.register_solv(
            solv_pdb=solv_pdb
        )
        solv_gro = os.path.basename(solv_pdb)+".pdb"
    soln.soln_setup(
        solute_itp_path=solute_itp,
        solute_gro_path=solute_gro,
        solute_itp_path_out=solute_itp_out,
        solv_gro_path=solv_gro,
        solv_itp_path=solv_itp,
        wf=wf,box_size=box_size
    )
    solute_name = os.path.basename(solute_itp_out)
    solute_name = os.path.splitext(solute_name)[0]
    solv_name = os.path.basename(solv_itp)
    solv_name = os.path.splitext(solv_name)[0]
    topfile = f"{solute_name}_{solv_name}_solvated.top"
    grofile = f"{solute_name}_{solv_name}_solvated.gro"
    soln_run = RunThrow(
        topfile=topfile,
        init_grofile=grofile
    )
    soln_run.run_EM()
    soln_run.run_nvt(T=T,run_ns=nvt_ns,grps=["MOL","SOL"])
    soln_run.run_npt(T=T,run_ns=npt_ns,grps=["MOL","SOL"])
    soln_run.run_production(T=T,run_ns=production_ns,grps=["MOL","SOL"])

