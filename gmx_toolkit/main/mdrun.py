import flet as ft

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
        elif value == "溶液のMD":
            self.page.clean()
        elif value == "WAXSプロファイル計算":
            self.page.clean()
        elif value == "密度計算":
            self.page.clean()

def main(page: ft.Page):
    GMXOperator(page)

ft.app(target=main)
