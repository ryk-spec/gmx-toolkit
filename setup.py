from setuptools import setup, find_packages

setup(
    name='gmx_toolkit',  # パッケージ名
    version='0.1.0',  # バージョン
    packages=find_packages(),  # パッケージを自動検出
    install_requires=[  # 必須の依存関係
        'numpy>=1.26.4',  # 例: numpyのバージョン指定
        'scipy>=1.12.1',
        'mdtraj>=1.9.1'
    ],
)
