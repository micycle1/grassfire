# tests2

This directory contains test code that was extracted from the main library files in `src/grassfire/`.

## Test Files

- `test_calc.py` - Tests for calc module
- `test_collapse.py` - Tests for collapse module (including main() and test functions)
- `test_line2d.py` - Tests for line2d module (including main() and helper test functions)
- `test_regression.py` - Tests for regression module
- `test_transform.py` - Tests for transform module
- `test_vectorops.py` - Tests for vectorops module

These tests were previously embedded in the library files as `if __name__ == "__main__"` blocks or test functions that were called from main().

## Running Tests

To run these tests, you'll need to have the grassfire package and its dependencies installed. Then you can run individual test files:

```bash
python -m grassfire.tests2.test_transform
```

Or use a test runner like pytest if available:

```bash
pytest src/grassfire/tests2/
```
