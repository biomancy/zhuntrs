import math
from pathlib import Path

import zhuntrs


def test_random_sequence_matches_c_reference() -> None:
    resources = Path(__file__).resolve().parent / "resources"
    sequence_path = resources / "random.txt"
    sequence = sequence_path.read_bytes().strip()

    reference = (resources / "random.txt.Z-SCORE").read_text(encoding="ascii")
    header, *rows = reference.splitlines()
    filename, length_text, mindn_text, maxdn_text = header.split()
    length, mindn, maxdn = int(length_text), int(mindn_text), int(maxdn_text)
    assert filename == sequence_path.name
    assert len(sequence) == length
    assert len(rows) == length

    predicted = zhuntrs.Engine(mindn, maxdn, wrap=True).predict(sequence)
    assert isinstance(predicted, zhuntrs.Prediction)
    assert len(predicted) == length
    starts, ends, scores, windows, conformations = (
        predicted.start,
        predicted.end,
        predicted.zhscore,
        predicted.sequence,
        predicted.conformation,
    )
    assert all(len(column) == length for column in (starts, ends, scores, windows, conformations))
    sequence_text = sequence.decode("ascii").upper()

    for index, row in enumerate(rows):
        _, _, expected_score, expected_conformation = row.split()
        expected_end = index + len(expected_conformation)
        assert starts[index] == index
        assert ends[index] == expected_end
        assert conformations[index] == expected_conformation
        assert math.isfinite(scores[index])
        assert abs(scores[index] - float(expected_score)) < 1e-12
        expected_window = "".join(
            sequence_text[position % length] for position in range(index, expected_end)
        )
        assert windows[index] == expected_window
