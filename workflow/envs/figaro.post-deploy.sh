conda activate $CONDA_PREFIX
pwd
cd ./workflow/envs/figaro
pwd
python3 setup.py bdist_wheel
pip3 install --force-reinstall dist/*.whl
