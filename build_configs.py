release32 = ["-DCMAKE_BUILD_TYPE=Release"]
debug32 = ["-DCMAKE_BUILD_TYPE=Debug"]
release32nolp = ["-DCMAKE_BUILD_TYPE=Release", "-DUSE_LP=NO"]
debug32nolp = ["-DCMAKE_BUILD_TYPE=Debug", "-DUSE_LP=NO"]
release64 = ["-DCMAKE_BUILD_TYPE=Release", "-DALLOW_64_BIT=True", "-DCMAKE_CXX_FLAGS='-m64'"]
debug64 = ["-DCMAKE_BUILD_TYPE=Debug",   "-DALLOW_64_BIT=True", "-DCMAKE_CXX_FLAGS='-m64'"]
release64nolp = ["-DCMAKE_BUILD_TYPE=Release", "-DALLOW_64_BIT=True", "-DCMAKE_CXX_FLAGS='-m64'", "-DUSE_LP=NO"]
debug64nolp = ["-DCMAKE_BUILD_TYPE=Debug",   "-DALLOW_64_BIT=True", "-DCMAKE_CXX_FLAGS='-m64'", "-DUSE_LP=NO"]
minimal = ["-DCMAKE_BUILD_TYPE=Release", "-DDISABLE_PLUGINS_BY_DEFAULT=YES"]


only_symbolic = ["-DCMAKE_BUILD_TYPE=Release", "-DDISABLE_PLUGINS_BY_DEFAULT=YES", "-DPLUGIN_SYMBOLIC_SEARCH_ENGINE_ENABLED=TRUE"]
release =  release64 + only_symbolic
debug = debug64 + only_symbolic

DEFAULT = "release"
DEBUG = "debug"
