import ctypes
import sys


class PyWorkerWrapper:
    def __init__(self, libName: str = "PythonWorker", libDirectory: str = "./lib"):
        # Get system name
        system = sys.platform
        # Check if windows
        if system == "win32":
            libPath = libDirectory + "/" + libName + ".dll"
        # Check if linux
        elif system == "linux":
            libPath = libDirectory + "/" + libName + ".so"
        # Check if mac
        elif system == "darwin":
            libPath = libDirectory + "/" + libName + ".dylib"
        else:
            print(f"Unsupported system: {system}")
            sys.exit(1)

        # Load library
        try:
            self.lib = ctypes.CDLL(libPath)
        except OSError as e:
            print(f"Error loading DLL: {e}")
            sys.exit(1)

        self.lib.PyWorker.argtypes = [
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
            ctypes.c_bool,  # writeInp
            ctypes.c_bool,  # writeMsh
            ctypes.c_bool,  # runFltk
            ctypes.c_bool,  # isPlaneStress
            ctypes.c_bool,  # verbose
        ]
        self.lib.PyWorker.restype = ctypes.c_double

        # Set default values
        self.L: float = 10
        self.B: float = 10
        self.ksi: float = 10
        self.a_b: float = 1
        self.loadValue: float = 1
        self.matrixModulus: float = 1
        self.matrixPoisson: float = 0.3
        self.inclusionModulus: float = 1
        self.inclusionPoisson: float = 0.2
        self.meshSize: float = 0.2
        self.refinementFactor: float = 1
        self.meshAlgorithm: int = 8
        self.isSerendipity = True
        self.convertToSquare = False
        self.writeInp = False
        self.writeMsh = False
        self.runFltk = False
        self.isPlaneStress = False
        self.verbose = False

    def run(self) -> float:
        return self.lib.PyWorker(
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
            self.writeInp,
            self.writeMsh,
            self.runFltk,
            self.isPlaneStress,
            self.verbose,
        )
