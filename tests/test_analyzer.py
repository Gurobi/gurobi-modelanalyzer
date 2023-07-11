from unittest import TestCase
import pathlib

import model_analyzer.analyzer as analyzer


class TestAnalyzer(TestCase):
    def setUp(self) -> None:
        here = pathlib.Path(__file__).parent
        self.file_path = str(here / "dataset" / "glass4.mps")

    def test_run(self):
        data = analyzer.run(self.file_path)
        print("Done")
