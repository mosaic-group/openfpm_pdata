source ~/openfpm_vars

CONDUIT_INSTALL_FOLDER=~/scratch/conduit
source ${CONDUIT_INSTALL_FOLDER}/.venv/bin/activate
export CONDUIT_DIR=${CONDUIT_INSTALL_FOLDER}/install-debug
export PYTHONPATH=${CONDUIT_INSTALL_FOLDER}/install-debug/python-modules/
