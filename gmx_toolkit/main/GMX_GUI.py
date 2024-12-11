import flet as ft
from gmx_toolkit.runner.gmx_run import solv_MD
class GMXOperator:
    def __init__(self, page: ft.Page):
        self.page = page
        self.page.title = "GUIでGROMACSの機能を使用できます"
        self.page.padding = 20
        
        # 初期化
        self.selected_value = ft.Ref[str]()
        self.init_ui()

    def init_ui(self):
        """初期UIのセットアップ"""
        # プルダウンメニュー
        self.dropdown = ft.Dropdown(
            ref=self.selected_value,
            label="処理を選んでください",
            options=[
                ft.dropdown.Option("溶媒のMD"),
                ft.dropdown.Option("溶液のMD"),
                ft.dropdown.Option("WAXSプロファイル計算"),
                ft.dropdown.Option("密度計算")
            ],
            width=300,
        )
        
        # 決定ボタン
        self.button = ft.ElevatedButton("決定", on_click=self.navigate)
        
        # ページのレイアウト
        self.page.add(
            ft.Column(
                controls=[self.dropdown,self.button],
                alignment=ft.MainAxisAlignment.CENTER,
                horizontal_alignment=ft.CrossAxisAlignment.CENTER,
            )
        )

    def navigate(self, e):
        """ページ遷移を処理"""
        value = self.selected_value.current.value
        if value == "溶媒のMD":
            self.page.clean()
            solvmd = SolvMD(page=self.page)
        elif value == "溶液のMD":
            self.page.clean()
        elif value == "WAXSプロファイル計算":
            self.page.clean()
        elif value == "密度計算":
            self.page.clean()

class SolvMD:

    def __init__(self,page: ft.Page):
        self.page = page
        self.solv_gro_picker = ft.FilePicker()
        self.solv_itp_picker = ft.FilePicker()
        self.solv_pdb_picker = ft.FilePicker()
        self.page.add(
            self.solv_gro_picker,
            self.solv_itp_picker,
            self.solv_pdb_picker
        )
        self.solv_gro_picker.on_result = self.solv_gro_picked
        self.solv_itp_picker.on_result = self.solv_itp_picked
        self.solv_pdb_picker.on_result = self.solv_pdb_picked
        self.solv_gro_field=ft.TextField(
                label="溶媒の.groファイル", width=500, height=50, multiline=True, min_lines=10,
                bgcolor=ft.colors.AMBER_100, border_radius=10)
        self.solv_gro_button = ft.ElevatedButton(
                "ファイルダイアログを開く", on_click=self.open_solv_gro,
                bgcolor=ft.colors.GREEN, color=ft.colors.AMBER_100, width=300, height=50,
        )

        self.solv_itp_field=ft.TextField(
                label="溶媒の.itpファイル", 
                width=500, height=50, multiline=True, min_lines=10,
                bgcolor=ft.colors.AMBER_100, border_radius=10)
        self.solv_itp_button = ft.ElevatedButton(
                "ファイルダイアログを開く", on_click=self.open_solv_itp,
                bgcolor=ft.colors.GREEN, color=ft.colors.AMBER_100, width=300, height=50,
        )
        

        self.solv_pdb_field=ft.TextField(
                label="溶媒の.pdbファイル(.groがないときに用いよ)", 
                width=500, height=50, multiline=True, min_lines=10,
                bgcolor=ft.colors.AMBER_100, border_radius=10)
        self.solv_pdb_button = ft.ElevatedButton(
                "ファイルダイアログを開く", on_click=self.open_solv_pdb,
                bgcolor=ft.colors.GREEN, color=ft.colors.AMBER_100, width=300, height=50,
        )
        self.box_size_field=ft.TextField(
                label="シミュレーションセルのサイズ in nm", width=500, height=50, multiline=True, min_lines=10,
                bgcolor=ft.colors.AMBER_100, border_radius=10)
        self.nvt_ns_field=ft.TextField(
                label="NVTランを何ns?", width=500, height=50, multiline=True, min_lines=10,
                bgcolor=ft.colors.AMBER_100, border_radius=10)
        self.npt_ns_field=ft.TextField(
                label="NPTランを何ns?", width=500, height=50, multiline=True, min_lines=10,
                bgcolor=ft.colors.AMBER_100, border_radius=10)
        self.production_ns_field=ft.TextField(
                label="Production ランを何ns?", width=500, height=50, multiline=True, min_lines=10,
                bgcolor=ft.colors.AMBER_100, border_radius=10)
        self.T_field=ft.TextField(
                label="温度は何K?", width=500, height=50, multiline=True, min_lines=10,
                bgcolor=ft.colors.AMBER_100, border_radius=10)
        self.run_button = ft.ElevatedButton(
                "MD開始", on_click=self.run,
                bgcolor=ft.colors.GREEN, color=ft.colors.AMBER_100, width=300, height=50,
        )
        self.nvt_ns_field.value = 0.5     
        self.npt_ns_field.value = 5
        self.production_ns_field.value = 1
        self.box_size_field.value = 5
        self.solv_pdb_field.value = "None"
        self.T_field.value = 301.15

        self.page.add(
            ft.Column(
                [
                    ft.Row([self.solv_gro_field,self.solv_gro_button]),
                    ft.Row([self.solv_itp_field,self.solv_itp_button]),
                    ft.Row([self.solv_pdb_field,self.solv_pdb_button]),
                    self.box_size_field,
                    self.nvt_ns_field,
                    self.npt_ns_field,
                    self.production_ns_field,
                    self.T_field,
                    self.run_button
                ]
            )
        )   

    def open_solv_gro(self, e):
        """ボタンを押した時に呼ばれる: ファイルピッカーを開く"""
        self.solv_gro_picker.pick_files()

    def solv_gro_picked(self, e):
        """ファイル選択が完了した時に呼ばれる"""
        self.solv_gro = self.solv_gro_picker.result.files[0].path
        self.solv_gro_field.value = self.solv_gro
        self.page.update()

    def open_solv_itp(self, e):
        """ボタンを押した時に呼ばれる: ファイルピッカーを開く"""
        self.solv_itp_picker.pick_files()

    def solv_itp_picked(self, e):
        """ファイル選択が完了した時に呼ばれる"""
        self.solv_itp = self.solv_itp_picker.result.files[0].path
        self.solv_itp_field.value = self.solv_itp
        self.page.update()

    def open_solv_pdb(self, e):
        """ボタンを押した時に呼ばれる: ファイルピッカーを開く"""
        self.solv_pdb_picker.pick_files()

    def solv_pdb_picked(self, e):
        """ファイル選択が完了した時に呼ばれる"""
        self.solv_pdb = self.solv_pdb_picker.result.files[0].path
        self.solv_pdb_field.value = self.solv_pdb
        self.page.update()

    def run(self,e):

        solv_MD(
            solv_gro=self.solv_gro_field.value,
            solv_itp=self.solv_itp_field.value,
            box_size=self.box_size_field.value,
            nvt_ns=self.nvt_ns_field.value,
            npt_ns=self.npt_ns_field.value,
            production_ns=self.production_ns_field.value,
            T=self.T_field.value,
            solv_pdb=self.solv_pdb_field.value
        )
        self.page.window.close()






def main(page: ft.Page):
    GMXOperator(page)

ft.app(target=main)
