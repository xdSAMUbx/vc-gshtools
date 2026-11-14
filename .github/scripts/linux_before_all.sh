#!/bin/sh
# Runs once inside the manylinux container before the wheels are built (cibuildwheel before_all).
# Installs Eigen (header-only, the only C++ dependency) and, when possible, ccache.
set -ex

# ccache is optional: if the package is not available, a pass-through shim keeps CMAKE_CXX_COMPILER_LAUNCHER=ccache working.
if ! yum install -y ccache; then
  printf "%s\n" "#!/bin/sh" "exec \"\$@\"" > /usr/local/bin/ccache
  chmod +x /usr/local/bin/ccache
  echo "ccache not available: using a pass-through shim"
fi

EIGEN_VERSION=3.4.0
cd /tmp
# Official GitLab archive of the 3.4.0 tag (retries: the download occasionally fails on shared runners).
curl -fsSL --retry 5 --retry-delay 5 -o eigen.tar.gz "https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz"
mkdir eigen-src
tar xzf eigen.tar.gz -C eigen-src --strip-components=1

/opt/python/cp312-cp312/bin/pip install cmake ninja
/opt/python/cp312-cp312/bin/cmake -S eigen-src -B eigen-build -G Ninja -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_TESTING=OFF
/opt/python/cp312-cp312/bin/cmake --install eigen-build
ls /usr/local/share/eigen3/cmake
