import ctypes

# 加载动态库
lib = ctypes.CDLL(
    "/Users/tengqing/Repository/myFEM/build/work/work.dylib"
)

# 函数参数类型设置
lib.work.argtypes = [
    ctypes.c_char_p,  # Path to the toml file
    ctypes.c_double,  # Elastic modulus of inclusion
    ctypes.c_double,  # ksi
    ctypes.c_double,  # a / b
    ctypes.c_bool,  # verbose
]

# 函数返回类型设置
lib.work.restype = ctypes.c_double

# 调用函数
tomlFilePath = "settings.toml".encode("utf-8")
result = lib.work(tomlFilePath, 1e-15, 10.0, 3.0, True)
print(f"Strain energy change: {result:.6f}")
