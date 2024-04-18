import libraryWrapper as lw
from Analytical import cylindricalInclusion

worker = lw.WorkerLibraryWrapper()

worker.B = 10
worker.ksi = 10
worker.a_b = 1
worker.matrixPoisson = 0.4
worker.inclusionModulus = 2
worker.refinementFactor = 5
worker.verbose = True
worker.runFltk = True
resultWorker = worker.run()
print("Strain energy change:")
print(f"  Integral on the interface result: {resultWorker:.6f}")
resultAnalytical = cylindricalInclusion(
    matrixPoisson=worker.matrixPoisson,
    inclusionModulus=worker.inclusionModulus,
    inclusionPoisson=worker.inclusionPoisson,
    inclusionCircleRadius=worker.B / worker.ksi,
)
print(f"  Analytical result: {resultAnalytical:.6f}")
print(
    f"Relative error: {abs(resultWorker - resultAnalytical) / abs(resultAnalytical)*100:.6f}%"
)
