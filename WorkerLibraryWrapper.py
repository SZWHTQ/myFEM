import ctypes


class WorkerLibraryWrapper:
    def __init__(self, path="./build/work/worker.dylib"):
        self.lib = ctypes.CDLL(path)
        # current folder via pwd
        self.lib.worker.argtypes = [
            ctypes.c_char_p,  # Path to the toml file
            ctypes.c_double,  # Elastic modulus of inclusion
            ctypes.c_double,  # ksi
            ctypes.c_double,  # a / b
            ctypes.c_bool,  # verbose
        ]
        self.lib.worker.restype = ctypes.c_double
        self.tomlFilePath = "work/settings.toml"
        self.inclusiveModulus = 1
        self.ksi = 10
        self.a_b = 1
        self.verbose = False

    def worker(self) -> float:
        return self.lib.worker(
            self.tomlFilePath.encode("utf-8"),
            self.inclusiveModulus,
            self.ksi,
            self.a_b,
            self.verbose,
        )
