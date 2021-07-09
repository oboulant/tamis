test:
	cd ./tests/c/ && $(MAKE) test && ./test
	python -m pytest

.PHONY: clean

clean:
	cd ./tests/c/ && $(MAKE) clean
