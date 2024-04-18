import libraryWrapper as lw

worker = lw.WorkerLibraryWrapper()

result = worker.run()
print(f"Strain energy change: {result:.6f}")
