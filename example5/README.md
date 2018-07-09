# GSoC 2018 MNE


### How to example 6:
replace the following files in your mne folder:

conf.py -> doc/conf.py
documentation.rst -> doc/documentation.rst
python_reference.rst -> doc/python_reference.rst
environment.yml -> environment.yml
environment2.yml -> environment2.yml
__init__.py -> mne/__init__.py
morph.py -> mne/morph.py

Afterwards navigate to the mne folder via terminal and execute:

python setup.py install

or

sudo -H python setup.py install

Then you should be able to use the newly introduced **SourceMorph** class
