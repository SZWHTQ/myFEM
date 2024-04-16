import ctypes

# 加载动态库
lib = ctypes.CDLL(
    "./build/work/worker.dylib"
)

# 函数参数类型设置
lib.worker.argtypes = [
    ctypes.c_char_p,  # Path to the toml file
    ctypes.c_double,  # Elastic modulus of inclusion
    ctypes.c_double,  # ksi
    ctypes.c_double,  # a / b
    ctypes.c_bool,  # verbose
]

# 函数返回类型设置
lib.worker.restype = ctypes.c_double

# 调用函数
tomlFilePath = "work/settings.toml".encode("utf-8")
result = lib.worker(tomlFilePath, 1e-15, 10.0, 1.0, False)
print(f"Strain energy change: {result:.6f}")
