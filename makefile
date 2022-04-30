develop :
	# export RUSTFLAGS='-L /opt/homebrew/include'
	maturin develop
	python test.py

release :
	# export RUSTFLAGS='-L /opt/homebrew/include'
	maturin build --release
	pip install /Users/lc663/GitHub/pyby/target/wheels/pyby-0.1.0-cp39-cp39-macosx_10_7_x86_64.whl --force-reinstall
	python test.py

