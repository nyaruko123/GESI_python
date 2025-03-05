import os
import sys


def StartupGCFB():
    """
    启动 GCFB 相关功能（类似于 MATLAB 中的 StartupGCFB）。
    可能需要加载 GCFB 相关的 Python 库或设置路径。
    """
    DirRoot = os.path.dirname(os.path.abspath(__file__))
    DirGCFB = os.path.join(DirRoot, "../../../GitHub_Public/gammachirp-filterbank/GCFBv234/")

    # 将 GCFB 路径添加到 sys.path
    if DirGCFB not in sys.path:
        sys.path.append(DirGCFB)

    print("GCFB 路径已添加:", DirGCFB)

    # 如果 GCFB 需要额外初始化，可以在这里调用相关代码
