"""Installed API checks; algorithm and input boundary cases live in Rust tests."""

import pickle
import re
import subprocess
import sys
import sysconfig
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor
from importlib.metadata import version
from pathlib import Path
from threading import Barrier
from typing import TypeVar, cast

import pytest

is_gil_enabled = cast(Callable[[], bool], getattr(sys, "_is_gil_enabled", lambda: True))
GIL_BEFORE_IMPORT = is_gil_enabled()
import zhuntrs

SEQUENCE = b"nATGCGCGCGGCATGCAT"
T = TypeVar("T")


def rows(predictions: Iterable[zhuntrs.Prediction]) -> list[tuple[int, int, float, str, str]]:
    return [
        row
        for prediction in predictions
        for row in zip(
            prediction.start,
            prediction.end,
            prediction.zhscore,
            prediction.sequence,
            prediction.conformation,
            strict=True,
        )
    ]


def parallel(function: Callable[[], T], workers: int = 4) -> list[T]:
    barrier = Barrier(workers)

    def run(_index: int) -> T:
        barrier.wait(timeout=30)
        return function()

    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(run, range(workers)))


def test_wheel_version_and_free_threaded_import() -> None:
    cargo = (Path(__file__).resolve().parents[1] / "Cargo.toml").read_text()
    match = re.search(r'^version = "([^"]+)"$', cargo, re.MULTILINE)
    assert match is not None
    assert version("zhuntrs") == match.group(1)
    if sysconfig.get_config_var("Py_GIL_DISABLED"):
        assert not GIL_BEFORE_IMPORT, "Run with the default GIL setting"
        assert not is_gil_enabled(), "Importing zhuntrs enabled the GIL"


def test_engine_configuration_is_read_only() -> None:
    engine = zhuntrs.Engine(1, 4)
    assert (engine.mindn, engine.maxdn, engine.threshold, engine.wrap) == (1, 4, 0.0, False)
    for name in ("mindn", "maxdn", "threshold", "wrap", "extra"):
        with pytest.raises(AttributeError):
            setattr(engine, name, 0)


def test_readme_example() -> None:
    prediction = zhuntrs.Engine(3, 6, 1000).predict(b"nATGCGCGCGGCATGC")
    assert isinstance(prediction, zhuntrs.Prediction)
    assert len(prediction) == 3
    assert prediction.start == (0, 1, 2)
    assert prediction.end == (10, 11, 10)
    assert [int(score) for score in prediction.zhscore] == [2220, 1057, 3428]
    assert prediction.sequence == ("NATGCGCGCG", "ATGCGCGCGG", "TGCGCGCG")
    assert prediction.conformation == ("ASASASASAS", "SASASASASA", "ASASASAS")
    assert (
        zhuntrs.Prediction(
            prediction.start,
            prediction.end,
            prediction.zhscore,
            prediction.sequence,
            prediction.conformation,
        )
        == prediction
    )
    with pytest.raises(AttributeError):
        prediction.start = ()  # type: ignore[misc]
    with pytest.raises(ValueError, match="equal lengths"):
        zhuntrs.Prediction([0], [], [], [], [])


def test_prediction_columns_are_cached_immutable_tuples() -> None:
    prediction = zhuntrs.Engine(1, 2).predict(SEQUENCE)
    for name in ("start", "end", "zhscore", "sequence", "conformation"):
        column: tuple[object, ...] = getattr(prediction, name)
        assert isinstance(column, tuple)
        assert getattr(prediction, name) is column
        with pytest.raises(TypeError):
            column[0] = None
    assert prediction != object()


def test_stream_returns_prediction_chunks() -> None:
    engine = zhuntrs.Engine(1, 4, 1000, wrap=True)
    iterator = engine.stream(SEQUENCE, chunk=3)

    assert isinstance(iterator, zhuntrs.PredictionIterator)
    assert iter(iterator) is iterator

    chunks = list(iterator)
    for prediction in chunks:
        assert isinstance(prediction, zhuntrs.Prediction)
        assert len(prediction) <= 3

    assert rows(chunks) == rows([engine.predict(SEQUENCE)])
    assert next(iterator, None) is None
    with pytest.raises(StopIteration):
        iterator.consume()


def test_engine_and_prediction_pickle_round_trip() -> None:
    engine = zhuntrs.Engine(2, 5, 1000, True)
    originals = (engine, engine.predict(SEQUENCE))
    script = """
import pickle
import sys
import zhuntrs

engine, prediction = pickle.loads(sys.stdin.buffer.read())
assert type(engine) is zhuntrs.Engine
assert type(prediction) is zhuntrs.Prediction
assert engine.predict(b"nATGCGCGCGGCATGCAT") == prediction
sys.stdout.buffer.write(pickle.dumps((engine, prediction)))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        input=pickle.dumps(originals),
        capture_output=True,
        check=True,
        timeout=30,
    )
    for restored, original in zip(pickle.loads(completed.stdout), originals, strict=True):
        assert type(restored) is type(original)
        assert restored == original


def test_concurrent_consume_is_atomic() -> None:
    engine = zhuntrs.Engine(2, 5, wrap=True)
    iterator = engine.stream(SEQUENCE, chunk=2)
    first = next(iterator)

    def consume() -> zhuntrs.Prediction | None:
        try:
            return iterator.consume()
        except StopIteration:
            return None

    remaining = [result for result in parallel(consume) if result is not None]
    assert len(remaining) == 1
    assert rows([first, *remaining]) == rows([engine.predict(SEQUENCE)])
