rm -R pyengnet.egg-info
rm -R dist
rm -R build
#python3 -m pip install --upgrade setuptools
#python3 setup.cfg sdist bdist_wheel
python setup.py install
#twine upload dist/*