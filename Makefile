# QMeku
PYTHON ?= python3

.PHONY: test clean
test:             ## run the tests (needs RDKit; the xtb / MOPAC tests are skipped if they are missing)
	$(PYTHON) -m unittest discover -s tests -v

clean:
	find qmeku tests benchmarks -type d -name __pycache__ -prune -exec rm -rf {} +
	rm -f reduce_info.log
