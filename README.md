# myFEM
For some homework

### Element Implementation
- Rectangle Serendipity 2D

### Depends on
1. [Gmsh](https://gmsh.info/)>=4.12.0
2. [tomlplusplus](https://github.com/marzer/tomlplusplus)
3. [Eigen](https://eigen.tuxfamily.org/)
4. [VTK](https://vtk.org/)

### Build
Configure the project with CMake
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
```
Build the project
```bash
cmake --build build -j
```
Or just build the pythonWorker library
```bash
cmake --build build -j --target pythonWorker
```
Go to the python directory
```bash
cd python
```
Before run the python script, you may need to modify the `libraryWrapper.py`. For example, if you build the library in the `python/library/build/` directory and decide to run python script in `python` directory, you should modify the `libraryWrapper.py` as follows
```python
class WorkerLibraryWrapper:
    def __init__(self, libDirectory: str = "./library/build/Release"):
        libDirectory = "./library/build"
        # Get system name
        ...
```
Now, run the python script
```bash
python Question2.py
```

### Result
#### Question 2
![This Project](image/Question2.png)
#### Single run
![This Project](image/StressXX.png)
![COMSOL](image/COMSOL_StressXX.png)
