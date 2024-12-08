# 溶媒追加手順

- `python frcmod_calc.py`で追加パラメータを算出
- 適宜名前を変え，`ffnonbonded.itp, ffbonded.itp, ffbonded2.itp`に追加
- `solv.itp`のatomsを編集．resnameは`SOL`，対応原子タイプを設定