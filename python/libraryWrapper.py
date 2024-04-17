import ctypes


class WorkerLibraryWrapper:
    def __init__(self, libPath="./worker.dylib"):
        self.lib = ctypes.CDLL(libPath)
        # current folder via pwd
        self.lib.worker.argtypes = [
            ctypes.c_double,  # L
            ctypes.c_double,  # B
            ctypes.c_double,  # ksi
            ctypes.c_double,  # a_b
            ctypes.c_double,  # loadValue
            ctypes.c_double,  # matrixModulus
            ctypes.c_double,  # matrixPoisson
            ctypes.c_double,  # inclusionModulus
            ctypes.c_double,  # inclusionPoisson
            ctypes.c_double,  # meshSize
            ctypes.c_double,  # refinementFactor
            ctypes.c_int,  # meshAlgorithm
            ctypes.c_bool,  # isSerendipity
            ctypes.c_bool,  # convertToSquare
            ctypes.c_bool,  # isPlaneStress
            ctypes.c_bool,  # verbose
        ]
        self.lib.worker.restype = ctypes.c_double
        self.L = 10
        self.B = 10
        self.ksi = 10
        self.a_b = 1
        self.loadValue = 1
        self.matrixModulus = 1
        self.matrixPoisson = 0.3
        self.inclusionModulus = 1
        self.inclusionPoisson = 0.2
        self.meshSize = 0.2
        self.refinementFactor = 1
        self.meshAlgorithm = 8
        self.isSerendipity = True
        self.convertToSquare = False
        self.isPlaneStress = False
        self.verbose = False

    def run(self) -> float:
        return self.lib.worker(
            self.L,
            self.B,
            self.ksi,
            self.a_b,
            self.loadValue,
            self.matrixModulus,
            self.matrixPoisson,
            self.inclusionModulus,
            self.inclusionPoisson,
            self.meshSize,
            self.refinementFactor,
            self.meshAlgorithm,
            self.isSerendipity,
            self.convertToSquare,
            self.isPlaneStress,
            self.verbose,
        )
