import libraryWrapper as lw

worker = lw.WorkerLibraryWrapper()

worker.runFltk = True
result = worker.run()
print(f"Strain energy change: {result:.6f}")
