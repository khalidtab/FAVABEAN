conda activate $CONDA_PREFIX
cd ./workflow/envs/figaro
python3 setup.py bdist_wheel
pip3 install --force-reinstall dist/*.whl
