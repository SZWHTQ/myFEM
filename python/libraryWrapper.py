import ctypes
import sys


class WorkerLibraryWrapper:
    def __init__(self, libDirectory: str = "./library/build/Release"):
        # Get system name
        system = sys.platform
        # Check if windows
        if system == "win32":
            libPath = libDirectory + "/PythonWorker.dll"
        # Check if linux
        elif system == "linux":
            libPath = libDirectory + "/PythonWorker.so"
        # Check if mac
        elif system == "darwin":
            libPath = libDirectory + "/PythonWorker.dylib"

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
            ctypes.c_bool,  # isPlaneStress
            ctypes.c_bool,  # verbose
        ]
        self.lib.PyWorker.restype = ctypes.c_double

        # Set default values
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
            self.isPlaneStress,
            self.verbose,
        )
