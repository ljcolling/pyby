develop :
	maturin develop
	python test.py

release :
	maturin build --release
	pip install ./target/wheels/pyby-0.1.0-cp39-cp39-macosx_10_7_x86_64.whl --force-reinstall
	python test.py

install :
	pip install ./target/wheels/pyby-0.1.0.tar.gz
	python test.py
